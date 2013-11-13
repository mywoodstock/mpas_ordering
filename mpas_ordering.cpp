/**
 *  Project:
 *
 *  File: mpas_ordering.cpp
 *  Created: Nov 12, 2013
 *  Modified: Wed 13 Nov 2013 11:44:34 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <algorithm>

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
	if(parts.empty()) {
		graph_partition_filename_ = "";
		num_partitions_ = 1;
	} else {
		graph_partition_filename_ = parts;
		num_partitions_ = 0;				// will be set later
	} // if-else

	if(!init()) {
		std::cerr << "error in init()" << std::endl;
		exit(1);
	} // if
} // MPASElementOrder::MPASElementOrder()


MPASElementOrder::~MPASElementOrder() { }


bool comp_umap_elements(const std::pair<int, MPASElementOrder::mpas_element_t>& a,
					const std::pair<int, MPASElementOrder::mpas_element_t>& b) {
	return (a.second < b.second);
} // comp_umap_elements()


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

	std::ifstream graph_f(graph_info_filename_.c_str());
	//int num_nodes, num_edges;
	//graph_f >> num_nodes; graph_f >> num_edges;
	std::string line;
	std::getline(graph_f, line);

	//if(num_nodes != num_cells_) {
	//	std::cerr << "error: number of cells do not match" << std::endl;
	//	return false;
	//} // if

	std::ifstream parts_f;
	if(num_partitions_ != 1) parts_f.open(graph_partition_filename_.c_str());

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
		ctemp = strtok((char*) line.c_str(), " ");
		int itemp;
		while(ctemp != NULL) {
			itemp = atoi(ctemp);
			cell.neighbor_list_.push_back(itemp - 1);
			ctemp = strtok(NULL, " ");
		} // while

		element_list_.push_back(cell);
	} // for

	graph_f.close();
	if(num_partitions_ != 1) parts_f.close();

	num_partitions_ = partition_list_.size();

	// sort the vector according to the current ordering (partition, ordering)
	std::sort(element_list_.begin(), element_list_.end());
	// generate the original index map
	generate_original_index_map();
	// reorder partition element lists according to current element ordering
	//reorder_partition_elements();
	// generate new order_index_ based on within partition index
	reindex_ordering_index();	// to convert order index to within partitions

	delete[] z_cells;
	delete[] y_cells;
	delete[] x_cells;

	return true;
} // MPASElementOrder::init()


//MPASElementOrder::set_mpas_element_t& MPASElementOrder::set_mpas_element_t::operator[](unsigned int i) {
//} // MPASElementOrder::set_mpas_element_t::operator[]()


// generate map from original index to index in element_list_, assume element_list_ is sorted
bool MPASElementOrder::generate_original_index_map() {
	original_index_map_.clear();
	unsigned int curr_index = 0;
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin(); ei != element_list_.end(); ++ ei)
		original_index_map_[(*ei).original_index_] = curr_index ++;
	return true;
} // MPASElementOrder::generate_original_index_map()


// reorder the element lists in partition_list_ according to the current element ordering
// ... is this really needed? probably not
//bool MPASElementOrder::reorder_partition_elements() {
//	return true;
//} // MPASElementOrder::reorder_partition_elements()


// reindex the ordering_index_ based on the order in partition_list_
bool MPASElementOrder::reindex_ordering_index() {
	for(partition_map_t::iterator pi = partition_list_.begin(); pi != partition_list_.end(); ++ pi) {
		unsigned int part_index = 0;
		for(vec_uint_t::iterator ei = (*pi).second.begin(); ei != (*pi).second.end(); ++ ei) {
			element_list_[original_index_map_[*ei]].ordering_index_ = part_index ++;
		} // for
	} // for
	return true;
} // MPASElementOrder::reindex_ordering_index()


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
			// ...
			return true;

		default:
			std::cerr << "error: unknown reordering specified" << std::endl;
			return false;
	} // switch

	return true;
} // MPASElementOrder::reorder_elements()


// write out the graph info and partition files out
bool MPASElementOrder::save_elements_order(std::string filename_prefix) {

	std::ofstream graph_f(filename_prefix.c_str());

	std::stringstream name;
	name << filename_prefix << ".part." << partition_list_.size();
	std::ofstream parts_f(name.str());

	graph_f << element_list_.size();
	unsigned int num_edges = 0;
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin();
			ei != element_list_.end(); ++ ei) {
		num_edges += (*ei).neighbor_list_.size();
		parts_f << (*ei).partition_num_ << std::endl;
	} // for
	num_edges /= 2;
	graph_f << "\t" << num_edges << std::endl;

	parts_f.close();

	for(vec_mpas_element_t::const_iterator ei = element_list_.begin();
			ei != element_list_.end(); ++ ei) {
		for(vec_uint_t::const_iterator ni = (*ei).neighbor_list_.begin();
				ni != (*ei).neighbor_list_.end(); ++ ni) {
			graph_f << "\t" << (original_index_map_[*ni] + 1);
		} // for
		graph_f << std::endl;
	} // for

	graph_f.close();

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


bool MPASElementOrder::print_all() {
	std::cout << "** num_cells: " << num_cells_ << std::endl;
	std::cout << "** num_partitions: " << num_partitions_ << std::endl;
	std::cout << "** element_list: " << std::endl;
	for(vec_mpas_element_t::iterator i = element_list_.begin(); i != element_list_.end(); ++ i) {
		std::cout << "    ";
		(*i).print_all();
	} // for
	std::cout << "** partition_list: " << std::endl;
	for(partition_map_t::iterator pi = partition_list_.begin(); pi != partition_list_.end(); ++ pi) {
		std::cout << "    " << (*pi).first << ":\t{ ";
		for(vec_uint_t::iterator ei = (*pi).second.begin(); ei != (*pi).second.end(); ++ ei) {
			std::cout << (*ei) << " ";
		} // for
		std::cout << "}" << std::endl;
	} // for

	return true;
} // MPASElementOrder::print_all()


bool MPASElementOrder::MPASElementData::print_all() const {
	std::cout << original_index_ << "\t" << partition_num_ << "\t" << ordering_index_ << "\t"
				<< "(" << x_coord_ << ", " << y_coord_ << ", " << z_coord_ << ")\t[ ";
	for(std::vector<unsigned int>::const_iterator i = neighbor_list_.cbegin();
			i != neighbor_list_.cend(); ++ i)
		std::cout << (*i) << " ";
	std::cout << "]" << std::endl;

	return true;
} // MPASElementOrder::MPASElementData::print_all()
