/**
 *  Project:
 *
 *  File: mpas_ordering_hilbert.cpp
 *  Created: Jan 16, 2014
 *  Modified: Thu 16 Jan 2014 03:58:27 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>

#include "mpas_ordering.hpp"

typedef MPASElementOrder::mpas_element_coord_pair_t mpas_element_coord_pair_t;


uint64_t MPASElementOrder::generate_hilbert_key(unsigned int x, unsigned int y, unsigned int z,
												unsigned int depth) {
	uint64_t key = 0;
	// get the last 'depth' bits, rest are discarded
	uint64_t x_64 = (uint64_t) x << (64 - depth);
	uint64_t y_64 = (uint64_t) y << (64 - depth);
	uint64_t z_64 = (uint64_t) z << (64 - depth);
	uint64_t mask = 1;	// leftmost bit is 1, rest are 0
	mask = mask << 63;

	// ignoring z ...

	int i = depth - 1;
//	uint64_t m2 = mask;
//	for(; i >= 0; -- i) {
//		if((x_64 & m2) || (y_64 & m2)) break;
//		x_64 = x_64 << 1;
//		y_64 = y_64 << 1;
//	} // for
	for(; i >= 0; -- i) {
		key = key << 1;
		if(mask - (x_64 & mask)) key = key | 1;

		key = key << 1;
		if((x_64 ^ y_64) & mask) key = key | 1;

		if(y_64 & mask) {
			if(mask - (x_64 & mask)) {
				x_64 = ~ x_64;
				y_64 = ~ y_64;
			}
			x_64 = x_64 ^ y_64;
			y_64 = x_64 ^ y_64;
			x_64 = x_64 ^ y_64;
		} // if

		x_64 = x_64 << 1;
		y_64 = y_64 << 1;
	} // for

	return key;
} // MPASElementOrder::generate_hilbert_key()


bool MPASElementOrder::reorder_elements_hilbert_sfc() {
//	if(on_sphere_) {
//		if(num_partitions_ < 8) {
//			// do some temporary partitions
//			// ...
//		} // if
		// make sure the element list is sorted on partitions
		std::sort(element_list_.begin(), element_list_.end());
		// generate partition pointers
		std::vector <unsigned int> part_begin;
		int prev_part = -1, curr_part = -1;
		for(unsigned int i = 0; i < num_cells_; ++ i) {
			curr_part = element_list_[i].partition_num_;
			if(curr_part != prev_part) part_begin.push_back(i);	// because they are sorted
			prev_part = curr_part;
		} // for
		part_begin.push_back(num_cells_);
//	} // if

	// get the minmax points for each partition
	std::vector <mpas_element_coord_pair_t> part_minmax;
	vec_mpas_element_t::const_iterator ei = element_list_.begin(), ei_begin, ei_end;
	for(int i = 0, j = 0; i < part_begin.size() - 1; ++ i) {
		ei_begin = ei;
		while(j < part_begin[i + 1]) { ++ j; ++ ei; }
		ei_end = ei;
		part_minmax.push_back(minmax(ei_begin, ei_end));
	} // for
	ei_begin = ei; ei_end = element_list_.end();
	part_minmax.push_back(minmax(ei_begin, ei_end));

	unsigned int new_index = 0;

	// for each partition, generate and assign hilbert order numbers to elements
	ei = element_list_.begin();
	for(int p = 0; p < part_begin.size() - 1; ++ p) {
		unsigned int j = part_begin[p]; ei_begin = ei;
		while(j < part_begin[p + 1]) { ++ j; ++ ei; }
		ei_end = ei;
		// project partition points onto a tangent plane first
		vec_mpas_element_coord_t projected_coords;
		project_points_to_plane(ei_begin, ei_end, part_minmax[p], projected_coords);

		std::map <uint64_t, unsigned int> hilbert_indices;	// keeps sorted on morton keys
		unsigned int depth = 21;		// max 21 bits in each of the 3 dims

		for(unsigned int e = part_begin[p]; e < part_begin[p + 1]; ++ e) {
			// generate morton number of element_list_[e]
			unsigned int index = e - part_begin[p];
			unsigned int x = (unsigned int) (projected_coords[index].x_);
			unsigned int y = (unsigned int) (projected_coords[index].y_);
			unsigned int z = (unsigned int) (projected_coords[index].z_);
			uint64_t key = generate_hilbert_key(x, y, z, depth);
			while(hilbert_indices.find(key) != hilbert_indices.end()) {
				std::cout << "warning: more than one points have the same morton key " << key
						<< ". this method fails here." << std::endl;
				++ key;
			} // if
			hilbert_indices[key] = e;
		} // for

		// generate the new indices for this partition elements
		for(std::map <uint64_t, unsigned int> :: iterator m = hilbert_indices.begin();
				m != hilbert_indices.end(); ++ m) {
			element_list_[(*m).second].ordering_index_ = new_index ++;
		} // for
	} // for

	// sort according to the hilbert numbers (ordering_index_)
	std::sort(element_list_.begin(), element_list_.end());

	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_morton_sfc_new()
