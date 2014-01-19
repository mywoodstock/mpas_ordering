/**
 *  Project:
 *
 *  File: mpas_ordering_peano.cpp
 *  Created: Jan 18, 2014
 *  Modified: Sat 18 Jan 2014 06:41:22 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>

#include "mpas_ordering.hpp"

typedef MPASElementOrder::mpas_element_coord_pair_t mpas_element_coord_pair_t;


uint64_t MPASElementOrder::generate_peano_key(unsigned int x, unsigned int y, unsigned int z,
												unsigned int depth = 13) {
	uint64_t key = 0;
	uint64_t x_3 = decimal_to_ternary(x);
	uint64_t y_3 = decimal_to_ternary(y);
	uint64_t z_3 = decimal_to_ternary(z);
	// get the last 'depth' digits
	uint64_t mask = ((uint64_t) 1 << 43) - 1;
	//uint64_t mask = (uint64_t) 8796093022207;
	x_3 = x_3 & mask;
	y_3 = y_3 & mask;
	z_3 = z_3 & mask;

	// here, depth is number of digits in ternary, not binary
	// ignoring z ...

	uint64_t div = 1000000000000;
	unsigned int xp = 0, yp = 0;
	for(int i = 12; i >= 0 && div > 0; -- i) {
		key *= 100;
		unsigned int x_digit = x_3 / div;
		unsigned int y_digit = y_3 / div;
		unsigned int t1 = x_digit, t2;
		yp += t1;
		if(yp % 2 == 0) t2 = y_digit;
		else t2 = 2 - y_digit;
		xp += t2;
		key += 10 * t1 + t2;
		x_3 = x_3 % div;
		y_3 = y_3 % div;
		div = div / 10;
	} // for

	return ternary_to_decimal(key);
} // MPASElementOrder::generate_peano_key()


bool MPASElementOrder::reorder_elements_peano_sfc() {
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

	// for each partition, generate and assign peano order numbers to elements
	ei = element_list_.begin();
	for(int p = 0; p < part_begin.size() - 1; ++ p) {
		unsigned int j = part_begin[p]; ei_begin = ei;
		while(j < part_begin[p + 1]) { ++ j; ++ ei; }
		ei_end = ei;
		// project partition points onto a tangent plane first
		vec_mpas_element_coord_t projected_coords;
		project_points_to_plane(ei_begin, ei_end, part_minmax[p], projected_coords);

		std::map <uint64_t, unsigned int> peano_indices;	// keeps sorted on morton keys
		unsigned int depth = 13;		// max 13 digits in each of the 3 dims

		for(unsigned int e = part_begin[p]; e < part_begin[p + 1]; ++ e) {
			// generate morton number of element_list_[e]
			unsigned int index = e - part_begin[p];
			unsigned int x = (unsigned int) (projected_coords[index].x_);
			unsigned int y = (unsigned int) (projected_coords[index].y_);
			unsigned int z = (unsigned int) (projected_coords[index].z_);
			uint64_t key = generate_peano_key(x, y, z, depth);
			//std::cout << p << "\t[ " << x << "\t" << y << "\t" << z << " ]\t" << key << std::endl;
			while(peano_indices.find(key) != peano_indices.end()) {
				std::cout << "warning: more than one points have the same morton key " << key
						<< ". this method fails here." << std::endl;
				++ key;
			} // if
			peano_indices[key] = e;
		} // for

		// generate the new indices for this partition elements
		for(std::map <uint64_t, unsigned int> :: iterator m = peano_indices.begin();
				m != peano_indices.end(); ++ m) {
			element_list_[(*m).second].ordering_index_ = new_index ++;
		} // for
	} // for

	// sort according to the peano numbers (ordering_index_)
	std::sort(element_list_.begin(), element_list_.end());

	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_morton_sfc_new()
