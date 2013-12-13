/**
 *  Project:
 *
 *  File: mpas_ordering_xyzsort.cpp
 *  Created: Dec 13, 2013
 *  Modified: Fri 13 Dec 2013 11:01:25 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>

#include "mpas_ordering.hpp"


struct XYZComp {
	bool operator()(const MPASElementOrder::MPASElementData& a,
					const MPASElementOrder::MPASElementData& b) {
		return (a.partition_num_ != b.partition_num_) ? (a.partition_num_ < b.partition_num_) :
				((a.x_coord_ != b.x_coord_) ? (a.x_coord_ < b.x_coord_) :
				((a.y_coord_ != b.y_coord_) ? (a.y_coord_ < b.y_coord_) :
				((a.z_coord_ != b.z_coord_) ? (a.z_coord_ < b.z_coord_) :
				false)));
	} // operator()
};


struct XComp {
	bool operator()(const MPASElementOrder::MPASElementData& a,
					const MPASElementOrder::MPASElementData& b) {
		return (a.partition_num_ != b.partition_num_) ? (a.partition_num_ < b.partition_num_) :
				(a.x_coord_ < b.x_coord_);
	} // operator()
};


struct YComp {
	bool operator()(const MPASElementOrder::MPASElementData& a,
					const MPASElementOrder::MPASElementData& b) {
		return (a.partition_num_ != b.partition_num_) ? (a.partition_num_ < b.partition_num_) :
				(a.y_coord_ < b.y_coord_);
	} // operator()
};


struct ZComp {
	bool operator()(const MPASElementOrder::MPASElementData& a,
					const MPASElementOrder::MPASElementData& b) {
		return (a.partition_num_ != b.partition_num_) ? (a.partition_num_ < b.partition_num_) :
				(a.z_coord_ < b.z_coord_);
	} // operator()
};


bool MPASElementOrder::reorder_elements_xyz_sort() {
	//ZComp z_comp;
	//std::sort(element_list_.begin(), element_list_.end(), z_comp);
	//YComp y_comp;
	//std::sort(element_list_.begin(), element_list_.end(), y_comp);
	//XComp x_comp;
	//std::sort(element_list_.begin(), element_list_.end(), x_comp);
	XYZComp xyz_comp;
	std::sort(element_list_.begin(), element_list_.end(), xyz_comp);
	return reindex_ordering_index();
} // MPASElementOrder::reorder_elements_xyz_sort()
