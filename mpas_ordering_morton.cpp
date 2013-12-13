/**
 *  Project:
 *
 *  File: mpas_ordering_morton.cpp
 *  Created: Dec 13, 2013
 *  Modified: Fri 13 Dec 2013 12:46:31 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>

#include "mpas_ordering.hpp"


// for a given range in element_list_, find min and max coordinates
std::pair <mpas_element_coord_t, mpas_element_coord_t>
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

	return std::pair <mpas_element_coord_t, mpas_element_coord_t> (min_point, max_point);
} // MPASElementOrder::minmax()

bool MPASElementOrder::reorder_elements_morton() {
	if(num_partitions_ == 1) {
		// do some temporary partitioning ...
	} // if
	// make sure the element list is sorted on partitions
	std::sort(element_list_.begin(), element_list_.end());
	// generate partition pointers
	std::vector <unsigned int> part_begin;
	int prev_part = curr_part = -1;
	for(unsigned int i = 0; i < num_cells_; ++ i) {
		curr_part = element_list_[i].partition_num_;
		if(curr_part != prev_part) part_begin.push_back(i);
		prev_part = curr_part;
	} // for

	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_xyz_sort()
