/**
 *  Project:
 *
 *  File: mpas_ordering.cpp
 *  Created: Nov 12, 2013
 *  Modified: Tue 12 Nov 2013 10:53:11 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>

#include "mpas_ordering.hpp"
#include "netcdf_utils.h"


MPASElementOrder::MPASElementOrder(std::string grid, std::string graph) {
	netcdf_grid_filename_ = grid;
	graph_info_filename_ = graph;
	graph_partition_filename_ = std::string("");
	num_partitions_ = 1;

	if(!init()) {
		std::cerr << "error in init()" << std::endl;
		exit(1);
	} // if
} // MPASElementOrder::MPASElementOrder()


MPASElementOrder::MPASElementOrder(std::string grid, std::string graph, std::string parts) {
	netcdf_grid_filename_ = grid;
	graph_info_filename_ = graph;
	graph_partition_filename_ = parts;
	num_partitions_ = 0;				// will be set later

	if(!init()) {
		std::cerr << "error in init()" << std::endl;
		exit(1);
	} // if
} // MPASElementOrder::MPASElementOrder()


bool MPASElementOrder::init() {
	// read in data from grid file
	// read in graph info
	// read in partitions
	
	num_cells_ = netcdf_mpas_read_dim(netcdf_grid_filename_, "nCells");

	double *x_cells, *y_cells, *z_cells;
	x_cells = new (std::nothrow) double[num_cells_];
	y_cells = new (std::nothrow) double[num_cells_];
	z_cells = new (std::nothrow) double[num_cells_];
	netcdf_mpas_read_xyzcell(netcdf_grid_filename_, num_cells_, x_cells, y_cells, z_cells);

	std::fstream graph_f(graph_info_filename_);
	int num_nodes, num_edges;
	graph_f >> num_nodes; graph_f >> num_edges;
	if(num_nodes != num_cells_) {
		std::cerr << "error: number of cells do not match" << std::endl;
		return false;
	} // if

	std::fstream parts_f;
	if(num_partitions_ != 1) parts_f.open(graph_partition_filename_);

	std::string line;
	for(int i = 0; i < num_cells_; ++ i) {
		MPASElementData cell;
		cell.original_index_ = i;
		cell.ordering_index_ = i;
		if(num_partitions_ != 1) parts_f >> cell.partition_num_;
		else cell.partition_num_ = 0;
		cell.x_coord_ = x_cells[i];
		cell.y_coord_ = y_cells[i];
		cell.z_coord_ = z_cells[i];
		partition_list_[cell.partition_num_].push_back(i);

		std::getline(graph_f, line);
		char* ctemp;
		ctemp = strtok(line.c_str(), " ");
		int itemp;
		while(ctemp != NULL) {
			itemp = atoi(ctemp);
			cell.neighbor_list_.push_back(itemp);
			ctemp = strtok(NULL, " ");
		} // while

		element_list_.push_back(cell);
	} // for

	graph_f.close();
	if(num_partitions_ != 1) parts_f.close();

	num_partitions_ = partition_list_.size();

	delete[] z_cells;
	delete[] y_cells;
	delete[] x_cells;

	return true;
} // MPASElementOrder::init()


bool MPASElementOrder::reorder_elements(sfc_t sfc) {
	switch(sfc) {
		case MORTON_SFC:
			return reorder_elements_morton_sfc();

		case HILBERT_SFC:
			return reorder_elements_hilbert_sfc();

		case XYZ_SORT:
			return reorder_elements_xyz_sort();

		case NONE:
			// back to original ordering
			for(std::vector<MPASElementData>::iterator i = element_list_.begin();
					i != element_list_.end(); ++ i) {
				(*i).ordering_index_ = (*i).original_index_;
			} // for
			return true;

		default:
			std::cerr << "error: unknown reordering specified" << std::endl;
			return false;
	} // switch

	return true;
} // MPASElementOrder::reorder_elements()


bool MPASElementOrder::save_elements_order(std::string filename_prefix) {

	return false;
} // MPASElementOrder::save_elements_order()


bool MPASElementOrder::reorder_elements_morton_sfc() {
	return false;
} // MPASElementOrder::reorder_elements_morton_sfc()


bool MPASElementOrder::reorder_elements_hilbert_sfc() {
	return false;
} // MPASElementOrder::reorder_elements_hilbert_sfc()


bool MPASElementOrder::reorder_elements_xyz_sort() {
	return false;
} // MPASElementOrder::reorder_elements_xyz_sort()

