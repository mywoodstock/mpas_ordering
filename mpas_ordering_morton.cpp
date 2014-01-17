/**
 *  Project:
 *
 *  File: mpas_ordering_morton.cpp
 *  Created: Dec 13, 2013
 *  Modified: Thu 16 Jan 2014 11:36:16 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>

#include "mpas_ordering.hpp"

typedef MPASElementOrder::mpas_element_coord_pair_t mpas_element_coord_pair_t;

// for a given range in element_list_, find min and max coordinates
/*mpas_element_coord_pair_t
MPASElementOrder::minmax(vec_mpas_element_t::const_iterator& begin,
						vec_mpas_element_t::const_iterator& end) const {
	double d_max = std::numeric_limits<double>::max();
	double d_min = std::numeric_limits<double>::min();
	mpas_element_coord_t min_point(d_max, d_max, d_max), max_point(d_min, d_min, d_min);
	for(vec_mpas_element_t::const_iterator i = begin; i != end; ++ i) {
		min_point.x_ = (min_point.x_ > (*i).x_coord_) ? (*i).x_coord_ : min_point.x_;
		max_point.x_ = (max_point.x_ < (*i).x_coord_) ? (*i).x_coord_ : max_point.x_;
		min_point.y_ = (min_point.y_ > (*i).y_coord_) ? (*i).y_coord_ : min_point.y_;
		max_point.y_ = (max_point.y_ < (*i).y_coord_) ? (*i).y_coord_ : max_point.y_;
		min_point.z_ = (min_point.z_ > (*i).z_coord_) ? (*i).z_coord_ : min_point.z_;
		max_point.z_ = (max_point.z_ < (*i).z_coord_) ? (*i).z_coord_ : max_point.z_;
	} // for
	return mpas_element_coord_pair_t(min_point, max_point);
} // MPASElementOrder::minmax()


// for a given range of coords, find min and max
mpas_element_coord_pair_t
MPASElementOrder::minmax(vec_mpas_element_coord_t::const_iterator& begin,
						vec_mpas_element_coord_t::const_iterator& end) const {
	vec_mpas_element_coord_t::const_iterator i = begin;
	mpas_element_coord_t min_point((*i).x_, (*i).y_, (*i).z_), max_point((*i).x_, (*i).y_, (*i).z_);
	++ i;
	for(; i != end; ++ i) {
		min_point.x_ = (min_point.x_ > (*i).x_) ? (*i).x_ : min_point.x_;
		max_point.x_ = (max_point.x_ < (*i).x_) ? (*i).x_ : max_point.x_;
		min_point.y_ = (min_point.y_ > (*i).y_) ? (*i).y_ : min_point.y_;
		max_point.y_ = (max_point.y_ < (*i).y_) ? (*i).y_ : max_point.y_;
		min_point.z_ = (min_point.z_ > (*i).z_) ? (*i).z_ : min_point.z_;
		max_point.z_ = (max_point.z_ < (*i).z_) ? (*i).z_ : max_point.z_;
	} // for
	return mpas_element_coord_pair_t(min_point, max_point);
} // MPASElementOrder::minmax()


bool pair_compare(std::pair <unsigned int, unsigned int> a, std::pair <unsigned int, unsigned int> b) {
	return a.first < b.first;
} // pair_compare


bool MPASElementOrder::project_points_to_plane(vec_mpas_element_t::const_iterator& begin,
							vec_mpas_element_t::const_iterator& end,
							mpas_element_coord_pair_t minmax_coords,
							vec_mpas_element_coord_t& projected_coords) {
	projected_coords.clear();
	mpas_element_coord_t center_point, projection;
	center_point.x_ = (minmax_coords.first.x_ + minmax_coords.second.x_) / 2;
	center_point.y_ = (minmax_coords.first.y_ + minmax_coords.second.y_) / 2;
	center_point.z_ = (minmax_coords.first.z_ + minmax_coords.second.z_) / 2;
	//std::cout << "min: (" << minmax_coords.first.x_ << ", "
	//					<< minmax_coords.first.y_ << ", "
	//					<< minmax_coords.first.z_ << ") max: ("
	//					<< minmax_coords.second.x_ << ", "
	//					<< minmax_coords.second.y_ << ", "
	//					<< minmax_coords.second.z_ << ") center: ("
	//					<< center_point.x_ << ", "
	//					<< center_point.y_ << ", "
	//					<< center_point.z_ << ")" << std::endl;
	double norm = center_point.x_ * center_point.x_ +
					center_point.y_ * center_point.y_ +
					center_point.z_ * center_point.z_;
	double a = center_point.x_ / norm;
	double b = center_point.y_ / norm;
	double c = center_point.z_ / norm;
	double abc = a * a + b * b + c * c;
	double k = c / sqrt(abc);
	double s = sqrt(1.0 - k * k);
	double l = 1.0 - k;
	double X = b / sqrt(a * a + b * b);
	double Y = - a / sqrt(a * a + b * b);
	double Z = 0.0;
	// rotation matrix
	double r00 = X * X * l + k,     r01 = X * Y * l - Z * s, r02 = X * Z * l + Y * s;
	double r10 = Y * X * l + Z * s, r11 = Y * Y * l + k,     r12 = Y * Z * l - X * s;
	double r20 = Z * X * l - Y * s, r21 = Z * Y * l + X * s, r22 = Z * Z * l + k;
	for(vec_mpas_element_t::const_iterator i = begin; i != end; ++ i) {
		double xp = (*i).x_coord_, yp = (*i).y_coord_, zp = (*i).z_coord_;
		double t = (1.0 - (a * xp + b * yp + c * zp)) / abc;
		double x1 = xp + a * t;
		double y1 = yp + b * t;
		double z1 = zp + c * t;
		projection.x_ = r00 * x1 + r01 * y1 + r02 * z1;
		projection.y_ = r10 * x1 + r11 * y1 + r12 * z1;
		projection.z_ = 0.0; //r20 * x1 + r21 * y1 + r22 * z1;
		projected_coords.push_back(projection);
		//std::cout << "*@@ [ " << xp << " , " << yp << " , " << zp << " ] -> [ "
		//			<< projection.x_ << " , " << projection.y_ << " , " << projection.z_ << " ]"
		//			<< std::endl;
	} // for

//	// make them all positive
//	vec_mpas_element_coord_t::const_iterator cb = projected_coords.begin();
//	vec_mpas_element_coord_t::const_iterator ce = projected_coords.end();
//	mpas_element_coord_pair_t new_minmax = minmax(cb, ce);
//	for(vec_mpas_element_coord_t::iterator i = projected_coords.begin(); i != projected_coords.end(); ++ i) {
//		if(new_minmax.first.x_ < 0.0) (*i).x_ = (*i).x_ - new_minmax.first.x_;
//		if(new_minmax.first.y_ < 0.0) (*i).y_ = (*i).y_ - new_minmax.first.y_;
//		//if(new_minmax.first.z_ < 0.0) (*i).z_ = (*i).z_ - new_minmax.first.z_;
//	} // for

	// sort on x and re-assign coords
	// sort on y and re-assign coords
	std::vector <std::pair <unsigned int, unsigned int> > x_mapper, y_mapper;
	unsigned int idx = 0;
	for(vec_mpas_element_coord_t::iterator i = projected_coords.begin();
		   	i != projected_coords.end(); ++ i, ++ idx) {
		x_mapper.push_back(std::pair <unsigned int, unsigned int> ((*i).x_, idx));
		y_mapper.push_back(std::pair <unsigned int, unsigned int> ((*i).y_, idx));
	} // for
	std::sort(x_mapper.begin(), x_mapper.end(), pair_compare);
	std::sort(y_mapper.begin(), y_mapper.end(), pair_compare);
	for(idx = 0; idx < projected_coords.size(); ++ idx) {
		projected_coords[x_mapper[idx].second].x_ = idx;
		projected_coords[y_mapper[idx].second].y_ = idx;
	} // for

	return true;
} // MPASElementOrder::project_points_to_plane()
*/

// given x, y and z coordinated as unsigned ints and height, generate the morton order key
uint64_t MPASElementOrder::generate_morton_key(unsigned int x, unsigned int y, unsigned int z,
												unsigned int height) {
	uint64_t key = 1;
	// get the last 'height' bits
	uint64_t x_64 = (uint64_t) x << (64 - height);
	uint64_t y_64 = (uint64_t) y << (64 - height);
	uint64_t z_64 = (uint64_t) z << (64 - height);
	uint64_t mask = 1;
	mask = mask << 63;	// leftmost bit is 1, rest are 0

	for(unsigned int i = 0; i < height; ++ i) {
		key = key << 1;
		if(x_64 & mask) key = key | 1;
		x_64 = x_64 << 1;

		key = key << 1;
		if(y_64 & mask) key = key | 1;
		y_64 = y_64 << 1;

		key = key << 1;
		if(z_64 & mask) key = key | 1;
		z_64 = z_64 << 1;
	} // for

	return key;
} // MPASElementOrder::generate_morton_key()


bool MPASElementOrder::reorder_elements_morton_sfc_new() {
//	if(on_sphere_) {
//		if(num_partitions_ < 8) {
			// do some temporary partitions
			// ...
//		} // if
		// make sure the element list is sorted on partitions
		std::sort(element_list_.begin(), element_list_.end());
		// generate partition pointers
		std::vector <unsigned int> part_begin;
		int prev_part = -1, curr_part = -1;
		for(unsigned int i = 0; i < num_cells_; ++ i) {
			curr_part = element_list_[i].partition_num_;
			if(curr_part != prev_part) part_begin.push_back(i);
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

	// for each partition, generate and assign morton numbers to elements
	ei = element_list_.begin();
	for(int p = 0; p < part_begin.size() - 1; ++ p) {
		unsigned int j = part_begin[p]; ei_begin = ei;
		while(j < part_begin[p + 1]) { ++ j; ++ ei; }
		ei_end = ei;
		vec_mpas_element_coord_t projected_coords;
		project_points_to_plane(ei_begin, ei_end, part_minmax[p], projected_coords);

		// number of bits required for each x, y and z = height of the tree = levels_ - 1
		// therefore, height of more than 21 (64/3) is not supported
		unsigned int height = 21;	// calculate this better using minmax coords ...
		std::map <uint64_t, unsigned int> morton_indices;	// keeps sorted on morton keys

		for(unsigned int e = part_begin[p]; e < part_begin[p + 1]; ++ e) {
			// generate morton number of element_list_[e]
			unsigned int index = e - part_begin[p];
			unsigned int x = (unsigned int) (projected_coords[index].x_);
			unsigned int y = (unsigned int) (projected_coords[index].y_);
			unsigned int z = (unsigned int) (projected_coords[index].z_);
			uint64_t key = generate_morton_key(x, y, z, height);
			//std::cout << "## " << p << "-" << e << ". ( " << x << " , " << y << " , " << z
			//			<< " ) -> " << key << std::endl;
			while(morton_indices.find(key) != morton_indices.end()) {
				std::cout << "warning: more than one points have the same morton key " << key
						<< ". this method fails here." << std::endl;
				++ key;
			} // if
			morton_indices[key] = e;
		} // for

		// generate the new indices for this partition elements
		for(std::map <uint64_t, unsigned int> :: iterator m = morton_indices.begin();
				m != morton_indices.end(); ++ m) {
			//std::cout << p << ". " << (*m).first << " -> " << (*m).second << std::endl;
			element_list_[(*m).second].ordering_index_ = new_index ++;
		} // for
	} // for

	// sort according to the morton numbers (ordering_index_)
	std::sort(element_list_.begin(), element_list_.end());

	//for(int i = 0; i < element_list_.size(); ++ i)
	//	std::cout << element_list_[i].ordering_index_ << "," << element_list_[i].original_index_ << "\t";
	//std::cout << std::endl;

	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_morton_sfc_new()


bool MPASElementOrder::reorder_elements_morton_sfc() {
//	if(on_sphere_ && num_partitions_ < 8) {
		// do some temporary partitioning ...
//	} // if
	// make sure the element list is sorted on partitions
	std::sort(element_list_.begin(), element_list_.end());
	// generate partition pointers
	std::vector <unsigned int> part_begin;
	int prev_part = -1, curr_part = -1;
	for(unsigned int i = 0; i < num_cells_; ++ i) {
		curr_part = element_list_[i].partition_num_;
		if(curr_part != prev_part) part_begin.push_back(i);
		prev_part = curr_part;
	} // for
	part_begin.push_back(num_cells_);

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

	// for each partition, generate and assign morton numbers to elements
	for(int p = 0; p < part_begin.size() - 1; ++ p) {
		// number of bits required for each x, y and z = height of the tree = levels_ - 1
		// therefore, height of more than 21 (64/3) is not supported
		unsigned int height = 21;	// calculate this better using minmax coords ...
		std::map <uint64_t, unsigned int> morton_indices;	// keeps sorted on morton keys
		mpas_element_coord_t min_point = part_minmax[p].first;
		mpas_element_coord_t translate(0.0, 0.0, 0.0);
		if(min_point.x_ < 0.0) translate.x_ = - min_point.x_;
		if(min_point.y_ < 0.0) translate.y_ = - min_point.y_;
		if(min_point.z_ < 0.0) translate.z_ = - min_point.z_;
		for(unsigned int e = part_begin[p]; e < part_begin[p + 1]; ++ e) {
			// generate morton number of element_list_[e]
			unsigned int x = (unsigned int) (element_list_[e].x_coord_ + translate.x_);
			unsigned int y = (unsigned int) (element_list_[e].y_coord_ + translate.y_);
			unsigned int z = (unsigned int) (element_list_[e].z_coord_ + translate.z_);
			uint64_t key = generate_morton_key(x, y, z, height);
			//std::cout << "(" << x << ", " << y << ", " << z << ") -> " << key << std::endl;
			if(morton_indices.find(key) != morton_indices.end()) {
				std::cerr << "error: more than one points have the same morton key " << key
						<< ". this method fails here." << std::endl;
				exit(1);
			} // if
			morton_indices[key] = e;
		} // for

		// generate the new indices for this partition elements
		for(std::map <uint64_t, unsigned int> :: iterator m = morton_indices.begin();
				m != morton_indices.end(); ++ m) {
			//std::cout << p << ". " << (*m).first << " -> " << (*m).second << std::endl;
			element_list_[(*m).second].ordering_index_ = new_index ++;
		} // for
	} // for

	// sort according to the morton numbers (ordering_index_)
	std::sort(element_list_.begin(), element_list_.end());

	//for(int i = 0; i < element_list_.size(); ++ i)
	//	std::cout << element_list_[i].ordering_index_ << "," << element_list_[i].original_index_ << "\t";
	//std::cout << std::endl;

	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_morton_sfc()
