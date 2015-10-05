/**
 *  Project:
 *
 *  File: mpas_ordering.cpp
 *  Created: Nov 12, 2013
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <algorithm>

#include <netcdfcpp.h>

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


bool MPASElementOrder::reorder_elements(sfc_t sfc) {
	// reorder cells
	switch(sfc) {
		case RANDOM:
			reorder_elements_random();
			break;

		case MORTON_SFC:
			reorder_elements_morton_sfc_new();
			break;

		case HILBERT_SFC:
			reorder_elements_hilbert_sfc();
			break;

		case PEANO_SFC:
			reorder_elements_peano_sfc();
			break;

		case XYZ_SORT:
			reorder_elements_xyz_sort();
			break;

		case NONE:
			// back to original ordering
			// ...
			return false;

		default:
			std::cerr << "error: unknown reordering specified" << std::endl;
			return false;
	} // switch
	generate_original_index_map();

	// reorder edges
	return reorder_edges();
} // MPASElementOrder::reorder_elements()


// write out the graph info and partition files out
bool MPASElementOrder::save_elements_order(std::string filename_prefix) {
	// create the mappers
	generate_original_index_map();

#ifdef DEBUG
	for(int i = 0; i < element_list_.size(); ++ i)
		std::cout << original_index_map_[i] + 1 << " -> " << i + 1 << std::endl;
#endif

	// reorder stuff in the grid and write a new grid file
  std::cout << "** reordering and writing new grid file ..." << std::endl;
	reorder_grid(filename_prefix);
  std::cout << "** writing new graph file ..." << std::endl;
	save_graph_info(filename_prefix);
  std::cout << "** writing new partition file ..." << std::endl;
	if(num_partitions_ > 1) save_partition_info(filename_prefix);

	return true;
} // MPASElementOrder::save_elements_order()


bool MPASElementOrder::save_graph_info(std::string prefix) {
	std::stringstream name;
	if(num_partitions_ > 1)
		name << prefix << "." << num_partitions_ << ".info";
	else
		name << prefix << ".info";
	std::ofstream graph_f(name.str());
	graph_f << element_list_.size() << "\t" << num_edges_ << std::endl;
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin();
			ei != element_list_.end(); ++ ei) {
		for(vec_uint_t::const_iterator ni = (*ei).neighbor_list_.begin();
				ni != (*ei).neighbor_list_.end(); ++ ni) {
			graph_f << "\t" << (current_index_map_[*ni] + 1);
		} // for
		graph_f << std::endl;
	} // for
	graph_f.close();

	return true;
} // MPASElementOrder::save_graph_info()


bool MPASElementOrder::save_partition_info(std::string prefix) {
	std::stringstream name;
	name << prefix << "." << num_partitions_ << ".info.part." << num_partitions_;
	std::ofstream parts_f(name.str());
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin();
			ei != element_list_.end(); ++ ei) {
		parts_f << (*ei).partition_num_ << std::endl;
	} // for
	parts_f.close();

	return true;
} // MPASElementOrder::save_partition_info()


bool MPASElementOrder::init() {
	// read in data from grid file
	// read in graph info
	// read in partitions
	
	num_cells_ = netcdf_mpas_read_dim(netcdf_grid_filename_, "nCells");
	num_edges_ = netcdf_mpas_read_dim(netcdf_grid_filename_, "nEdges");
#ifdef DEBUG
	std::cout << "NUM CELLS = " << num_cells_ << ", NUM EDGES = " << num_edges_ << std::endl;
#endif

	// read cells spatial data
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

	//unsigned int num_edges = 0;
	for(int i = 0; i < num_cells_; ++ i) {
		mpas_element_t cell;
		cell.original_index_ = i;
		cell.ordering_index_ = i;
		if(num_partitions_ != 1) parts_f >> cell.partition_num_;
		else cell.partition_num_ = 0;
		cell.x_coord_ = x_cells[i];
		cell.y_coord_ = y_cells[i];
		cell.z_coord_ = z_cells[i];
		//partition_list_[cell.partition_num_].push_back(i);

		std::getline(graph_f, line);
		char* ctemp;
		ctemp = strtok((char*) line.c_str(), " ");
		int itemp;
		while(ctemp != NULL) {
			itemp = atoi(ctemp);
			cell.neighbor_list_.push_back(itemp - 1);
			ctemp = strtok(NULL, " ");
		} // while
		//num_edges += cell.neighbor_list_.size();
		element_list_.push_back(cell);
	} // for
	//num_edges = num_edges >> 1;

	delete[] z_cells;
	delete[] y_cells;
	delete[] x_cells;

	//if(num_edges != num_edges_) {
	//	std::cerr << "error: number of edges do not match: " << num_edges << " != " << num_edges_
	//				<< std::endl;
	//	return false;
	//} // if

	graph_f.close();
	if(num_partitions_ != 1) parts_f.close();

	//num_partitions_ = partition_list_.size();

	// sort the vector according to the current ordering (partition, ordering)
	std::sort(element_list_.begin(), element_list_.end());
	// generate the original index map
	generate_original_index_map();
	// reorder partition element lists according to current element ordering
	//reorder_partition_elements();
	// generate new order_index_ based on within partition index
	//reindex_ordering_index();	// to convert order index to within partitions
	generate_partition_list();

	// read edge data (for reordering edges later)
	int *cellsonedge = new (std::nothrow) int[num_edges_ * 2];	// 2 cells per edge
	netcdf_mpas_read_cellsonedge(netcdf_grid_filename_, num_edges_, cellsonedge);

	for(int i = 0; i < num_edges_; ++ i) {
		mpas_edge_t edge;
		edge.original_index_ = i;
		edge.ordering_index_ = i;
		edge.cell_index_[0] = cellsonedge[2 * i];
		edge.cell_index_[1] = cellsonedge[2 * i + 1];
		if(edge.cell_index_[0] != 0 && edge.cell_index_[1] != 0) {
			edge.partition_num_ = std::min(
							element_list_[current_index_map_[edge.cell_index_[0] - 1]].partition_num_,
							element_list_[current_index_map_[edge.cell_index_[1] - 1]].partition_num_);
		} else {
			if(edge.cell_index_[0] != 0) {
				edge.partition_num_ =
							element_list_[current_index_map_[edge.cell_index_[0] - 1]].partition_num_;
			} else {
				edge.partition_num_ =
							element_list_[current_index_map_[edge.cell_index_[1] - 1]].partition_num_;
			} // if-else
		} // if-else
		edge_list_.push_back(edge);
	} // for

	delete[] cellsonedge;

#ifdef DEBUG
	print_all();
#endif

	return true;
} // MPASElementOrder::init()


// generate maps from/to current index in element_list_ to/from original index
bool MPASElementOrder::generate_original_index_map() {
	original_index_map_.clear();
	current_index_map_.clear();
	unsigned int curr_index = 0;
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin(); ei != element_list_.end(); ++ ei) {
		original_index_map_[curr_index] = (*ei).original_index_;
		current_index_map_[(*ei).original_index_] = curr_index;
		++ curr_index;
	} // for
	return true;
} // MPASElementOrder::generate_original_index_map()


bool MPASElementOrder::generate_partition_list() {
	partition_list_.clear();
	for(vec_mpas_element_t::const_iterator ei = element_list_.begin(); ei != element_list_.end(); ++ ei) {
		partition_list_[(*ei).partition_num_].push_back((*ei).original_index_);
	} // for
	num_partitions_ = partition_list_.size();
	return true;
} // MPASElement::generate_partition_list()


// reindex the ordering_index_ within each partition according to current order
bool MPASElementOrder::reindex_ordering_index() {
	vec_uint_t pcounters(num_partitions_, 0);
	for(vec_mpas_element_t::iterator ei = element_list_.begin(); ei != element_list_.end(); ++ ei) {
		(*ei).ordering_index_ = pcounters[(*ei).partition_num_] ++;
	} // for
	return true;
} // MPASElementOrder::reindex_ordering_index()


// reindex and reorder edge list according to current cell ordering
bool MPASElementOrder::reorder_edges() {
	// assign new indices
	for(vec_mpas_edge_t::iterator ei = edge_list_.begin(); ei != edge_list_.end(); ++ ei) {
		if((*ei).cell_index_[0] != 0 && (*ei).cell_index_[1] != 0) {
			mpas_element_t cell0 = element_list_[current_index_map_[(*ei).cell_index_[0] - 1]];
			mpas_element_t cell1 = element_list_[current_index_map_[(*ei).cell_index_[1] - 1]];
			(*ei).partition_num_ = std::min(cell0.partition_num_, cell1.partition_num_);
			if(cell0.partition_num_ == cell1.partition_num_) {
				(*ei).partition_num_ = cell0.partition_num_;
				(*ei).ordering_index_ = std::min(cell0.ordering_index_, cell1.ordering_index_);
			} else if (cell0.partition_num_ < cell1.partition_num_) {
				(*ei).partition_num_ = cell0.partition_num_;
				(*ei).ordering_index_ = cell0.ordering_index_;
			} else {
				(*ei).partition_num_ = cell1.partition_num_;
				(*ei).ordering_index_ = cell1.ordering_index_;
			} // if-else
		} else if((*ei).cell_index_[0] != 0) {
			mpas_element_t cell0 = element_list_[current_index_map_[(*ei).cell_index_[0] - 1]];
			(*ei).partition_num_ = cell0.partition_num_;
			(*ei).ordering_index_ = cell0.ordering_index_;
		} else {
			mpas_element_t cell1 = element_list_[current_index_map_[(*ei).cell_index_[1] - 1]];
			(*ei).partition_num_ = cell1.partition_num_;
			(*ei).ordering_index_ = cell1.ordering_index_;
		} // if-else
	} // for

	// sort
	std::sort(edge_list_.begin(), edge_list_.end());

	// reindex the new indices
	vec_mpas_edge_t::iterator ei = edge_list_.begin();
	unsigned int prev_p = (*ei).partition_num_, prev_i = (*ei).ordering_index_, counter = 0;
	(*ei).ordering_index_ = counter ++;
	++ ei;
	for(; ei != edge_list_.end(); ++ ei) {
		if(prev_p != (*ei).partition_num_) counter = 0;
		(*ei).ordering_index_ = counter ++;
		prev_p = (*ei).partition_num_;
	} // for

	return generate_edge_original_index_map();
} // MPASElemetOrder::reorder_edges()


// generate maps from/to current index in edge_list_ to/from original index
bool MPASElementOrder::generate_edge_original_index_map() {
	edge_original_index_map_.clear();
	edge_current_index_map_.clear();
	unsigned int curr_index = 0;
	for(vec_mpas_edge_t::const_iterator ei = edge_list_.begin(); ei != edge_list_.end(); ++ ei) {
		edge_original_index_map_[curr_index] = (*ei).original_index_;
		edge_current_index_map_[(*ei).original_index_] = curr_index;
		++ curr_index;
	} // for
#ifdef DEBUG
	print_all();
#endif
	return true;
} // MPASElementOrder::generate_original_index_map()


bool MPASElementOrder::reorder_grid(std::string prefix) {
	std::stringstream grid_name;
	if(num_partitions_ > 1)
		grid_name << prefix << "." << num_partitions_ << ".nc";
	else
		grid_name << prefix << ".nc";
	#ifdef _64BITOFFSET
		NcFile ncid(netcdf_grid_filename_.c_str(), NcFile::ReadOnly, NULL, 0, NcFile::Offset64Bits);
		NcFile out_ncid(grid_name.str().c_str(), NcFile::Replace, NULL, 0, NcFile::Offset64Bits);
	#else
		NcFile ncid(netcdf_grid_filename_.c_str(), NcFile::ReadOnly);
		NcFile out_ncid(grid_name.str().c_str(), NcFile::Replace);
	#endif

	// read dims
	for(int i = 0; i < ncid.num_dims(); ++ i) {
		NcDim *dim_id = ncid.get_dim(i);
		// check if unlimited
		if((*dim_id).is_unlimited())
			out_ncid.add_dim((*dim_id).name());
		else
			out_ncid.add_dim((*dim_id).name(), (*dim_id).size());
	} // for

	// read in attributes
	for(int i = 0; i < ncid.num_atts(); ++ i) {
		NcAtt *att_id = ncid.get_att(i);
		// NOTE: assuming that all non-string attributes have only one value
		ncbyte att_byte;
		char *att_chars;
		short att_short;
		int att_int;
		float att_float;
		double att_double;
		int j = 0;
		switch((*att_id).type()) {
			case ncChar:
				att_chars = new (std::nothrow) char[(*att_id).num_vals()];
				for(; j < (*att_id).num_vals(); ++ j) att_chars[j] = (*att_id).as_char(j);
				att_chars[j] = '\0';
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), att_chars);
				break;

			case ncShort:
				att_short = (*att_id).as_short(0);
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), &att_short);
				break;

			case ncInt:
				att_int = (*att_id).as_int(0);
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), &att_int);
				break;

			case ncFloat:
				att_float = (*att_id).as_float(0);
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), &att_float);
				break;

			case ncDouble:
				att_double = (*att_id).as_double(0);
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), &att_double);
				break;

			case ncByte:
				att_byte = (*att_id).as_ncbyte(0);
				out_ncid.add_att((*att_id).name(), (*att_id).num_vals(), &att_byte);
				break;

			default:
				std::cerr << "error: unknown type for attribute" << std::endl;
				exit(1);
		} // switch
	} // for

	NcToken cell_dim_token = ncid.get_dim("nCells")->name();
	NcToken cell_id_token = ncid.get_var("indexToCellID")->name();
	NcToken cells_on_edge_token = ncid.get_var("cellsOnEdge")->name();
	NcToken cells_on_cell_token = ncid.get_var("cellsOnCell")->name();
	NcToken cells_on_vertex_token = ncid.get_var("cellsOnVertex")->name();

	NcToken edge_dim_token = ncid.get_dim("nEdges")->name();
	NcToken edge_id_token = ncid.get_var("indexToEdgeID")->name();
	NcToken edges_on_cell_token = ncid.get_var("edgesOnCell")->name();
	NcToken edges_on_edge_token = ncid.get_var("edgesOnEdge")->name();
	NcToken edges_on_vertex_token = ncid.get_var("edgesOnVertex")->name();

	// read vars and reorder data
	NcDim *dims[5];
	long dim_sizes[5];
	for(int i = 0; i < ncid.num_vars(); ++ i) {
		NcVar *var_id = ncid.get_var(i);
		int num_dims = (*var_id).num_dims();
		std::map <NcToken, int> var_dims;
		bool cells_to_reorder = false;
		int cell_dim_num = -1;
		bool edges_to_reorder = false;
		int edge_dim_num = -1;
		int size = 1;
		for(int j = 0; j < num_dims; ++ j) {
			dims[j] = (*var_id).get_dim(j);
			dim_sizes[j] = (*dims[j]).size();
			size *= (*dims[j]).size();
			var_dims[(*dims[j]).name()] = (*dims[j]).size();
			if((*dims[j]).name() == cell_dim_token && (*var_id).name() != cell_id_token) {
				cells_to_reorder = true;
				cell_dim_num = j;
			} // if
			if((*dims[j]).name() == edge_dim_token && (*var_id).name() != edge_id_token) {
				edges_to_reorder = true;
				edge_dim_num = j;
			} // if	
		} // for
		
#ifdef DEBUG
		std::cout << "*** Variable: " << (*var_id).name();
		std::cout << ", Cell reordering: " << (cells_to_reorder ? "yes" : "no");
		std::cout << ", Edge reordering: " << (edges_to_reorder ? "yes" : "no");
		std::cout << std::endl;
#endif

		// NOTE: assuming that at most one dimension may be equal to nCells
		// NOTE: assuming that no variable has more than 5 dimensions
		NcVar *out_var_id = out_ncid.add_var((*var_id).name(), (*var_id).type(),
												num_dims, (const NcDim**)dims);

		int *int_buff_in = NULL, *int_buff_out = NULL;
		double *dbl_buff_in = NULL, *dbl_buff_out = NULL;
		char *chr_buff = NULL;
		switch((*var_id).type()) {
			case ncByte:
			case ncShort:
			case ncFloat:
				std::cerr << "error: unimplemented variable type" << std::endl;
				exit(1);

			case ncChar:
				chr_buff = new (std::nothrow) char[size];
				(*var_id).get(chr_buff, dim_sizes);
				(*out_var_id).put(chr_buff, dim_sizes);
				delete[] chr_buff;
				break;

			case ncInt:
				int_buff_in = new (std::nothrow) int[size];
				(*var_id).get(int_buff_in, dim_sizes);
				if(cells_to_reorder) {	// cells
					int_buff_out = new (std::nothrow) int[size];
					reorder_data<int>(int_buff_in, int_buff_out, cell_dim_num, num_dims, dim_sizes);
					if(edges_to_reorder) {	// both
						int* temp = int_buff_out;
						int_buff_out = int_buff_in;
						int_buff_in = temp;
						reorder_edge_data<int>(int_buff_in, int_buff_out, edge_dim_num, num_dims, dim_sizes);
					} // if
				} else {
					if(edges_to_reorder) {	// only edges
						int_buff_out = new (std::nothrow) int[size];
						reorder_edge_data<int>(int_buff_in, int_buff_out, edge_dim_num, num_dims, dim_sizes);
					} else {	// none
						int_buff_out = int_buff_in;
					} // if-else
				} // if-else
				// update the cell numbers in cellsOnEdge, cellsOnCell, cellsOnVertex
				if((*var_id).name() == cells_on_edge_token ||
						(*var_id).name() == cells_on_cell_token ||
						(*var_id).name() == cells_on_vertex_token) {
					for(unsigned int i = 0; i < size; ++ i)
						if(int_buff_out[i] > 0)
							int_buff_out[i] = current_index_map_[int_buff_out[i] - 1] + 1;
				} // if
				if((*var_id).name() == edges_on_cell_token ||
						(*var_id).name() == edges_on_edge_token ||
						(*var_id).name() == edges_on_vertex_token) {
					for(unsigned int i = 0; i < size; ++ i)
						if(int_buff_out[i] > 0)
							int_buff_out[i] = edge_current_index_map_[int_buff_out[i] - 1] + 1;
				} // if
				(*out_var_id).put(int_buff_out, dim_sizes);
				if(cells_to_reorder || edges_to_reorder) delete[] int_buff_out;
				delete[] int_buff_in;
				break;

			case ncDouble:
				dbl_buff_in = new (std::nothrow) double[size];
				(*var_id).get(dbl_buff_in, dim_sizes);
				if(cells_to_reorder) {
					dbl_buff_out = new (std::nothrow) double[size];
					reorder_data<double>(dbl_buff_in, dbl_buff_out, cell_dim_num, num_dims, dim_sizes);
					(*out_var_id).put(dbl_buff_out, dim_sizes);
					delete[] dbl_buff_out;
				} else {
					if(edges_to_reorder) {
						dbl_buff_out = new (std::nothrow) double[size];
						reorder_edge_data<double>(dbl_buff_in, dbl_buff_out, edge_dim_num, num_dims, dim_sizes);
						(*out_var_id).put(dbl_buff_out, dim_sizes);
						delete[] dbl_buff_out;
					} else {
						(*out_var_id).put(dbl_buff_in, dim_sizes);
					} // if-else
				} // if-else
				delete[] dbl_buff_in;
				break;

			default:
				std::cerr << "error: unknown variable type" << std::endl;
				exit(1);
		} // switch

	} // for
	
	out_ncid.close();
	ncid.close();

	return true;
} // MPASElementOrder::reorder_grid()


// to reorder cell data
template<typename T>
bool MPASElementOrder::reorder_data(T* in, T* out, int cell_dim_num, int num_dims, long* dim_sizes) {
	// NOTE: assuming no more than three dimensions
	// NOTE: assuming nCells is either first or second dimension
	unsigned int offset = 0;
	switch(num_dims) {
		case 1:
			for(unsigned int i = 0; i < num_cells_; ++ i) {
				out[i] = in[original_index_map_[i]];
			} // for
			break;

		case 2:
			switch(cell_dim_num) {
				case 0:
					offset = dim_sizes[1];
					for(unsigned int i = 0; i < num_cells_; ++ i) {
						memcpy(&out[i * offset], &in[original_index_map_[i] * offset],
								offset * sizeof(T));
					} // for
					break;

				case 1:
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							out[i * dim_sizes[1] + j] = in[i * dim_sizes[1] + original_index_map_[j]];
						} // for
					} // for
					break;

				default:
					std::cerr << "error: mismatch in num_dims and cell_dim_num" << std::endl;
					exit(1);
			} // switch
			break;

		case 3:
			switch(cell_dim_num) {
				case 0:
					offset = dim_sizes[1] * dim_sizes[2];
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						memcpy(&out[i * offset], &in[original_index_map_[i] * offset],
								offset * sizeof(T));
					} // for
					break;

				case 1:
					offset = dim_sizes[2];
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							memcpy(&out[(i * dim_sizes[1] + j) * offset],
									&in[(i * dim_sizes[1] + original_index_map_[j]) * offset],
									offset * sizeof(T));
						} // for
					} // for
					break;

				case 2:
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							offset = (i * dim_sizes[1] + j) * dim_sizes[2];
							for(unsigned int k = 0; k < dim_sizes[2]; ++ k) {
								out[offset + k] = in[offset + original_index_map_[k]];
							} // for
						} // for
					} // for
					break;

				default:
					std::cerr << "error: mismatch in num_dims and cell_dim_num" << std::endl;
					exit(1);
			} // switch
			break;

		default:
			std::cerr << "error: case with more than 3 dimensions not implemented" << std::endl;
			exit(1);
	} // switch
	return true;
} // MPASElementOrder::reorder_data()


// to reorder edge data
template<typename T>
bool MPASElementOrder::reorder_edge_data(T* in, T* out, int edge_dim_num, int num_dims, long* dim_sizes) {
	// NOTE: assuming no more than three dimensions
	// NOTE: assuming nCells is either first or second dimension
	unsigned int offset = 0;
	switch(num_dims) {
		case 1:
			for(unsigned int i = 0; i < num_edges_; ++ i) {
				out[i] = in[edge_original_index_map_[i]];
			} // for
			break;

		case 2:
			switch(edge_dim_num) {
				case 0:
					offset = dim_sizes[1];
					for(unsigned int i = 0; i < num_edges_; ++ i) {
						memcpy(&out[i * offset], &in[edge_original_index_map_[i] * offset],
								offset * sizeof(T));
					} // for
					break;

				case 1:
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							out[i * dim_sizes[1] + j] = in[i * dim_sizes[1] + edge_original_index_map_[j]];
						} // for
					} // for
					break;

				default:
					std::cerr << "error: mismatch in num_dims and edge_dim_num" << std::endl;
					exit(1);
			} // switch
			break;

		case 3:
			switch(edge_dim_num) {
				case 0:
					offset = dim_sizes[1] * dim_sizes[2];
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						memcpy(&out[i * offset], &in[edge_original_index_map_[i] * offset],
								offset * sizeof(T));
					} // for
					break;

				case 1:
					offset = dim_sizes[2];
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							memcpy(&out[(i * dim_sizes[1] + j) * offset],
									&in[(i * dim_sizes[1] + edge_original_index_map_[j]) * offset],
									offset * sizeof(T));
						} // for
					} // for
					break;

				case 2:
					for(unsigned int i = 0; i < dim_sizes[0]; ++ i) {
						for(unsigned int j = 0; j < dim_sizes[1]; ++ j) {
							offset = (i * dim_sizes[1] + j) * dim_sizes[2];
							for(unsigned int k = 0; k < dim_sizes[2]; ++ k) {
								out[offset + k] = in[offset + edge_original_index_map_[k]];
							} // for
						} // for
					} // for
					break;

				default:
					std::cerr << "error: mismatch in num_dims and edge_dim_num" << std::endl;
					exit(1);
			} // switch
			break;

		default:
			std::cerr << "error: case with more than 3 dimensions not implemented" << std::endl;
			exit(1);
	} // switch
	return true;
} // MPASElementOrder::reorder_data()


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
	std::cout << "** edge_list:" << std::endl;
	for(vec_mpas_edge_t::const_iterator i = edge_list_.begin(); i != edge_list_.end(); ++ i) {
		std::cout << "    ";
		(*i).print_all();
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


bool MPASElementOrder::MPASEdgeData::print_all() const {
	std::cout << original_index_ << "\t" << partition_num_ << "\t" << ordering_index_ << "\t"
				<< "[" << cell_index_[0] << ", " << cell_index_[1] << "]" << std::endl;
	return true;
} // MPASElementOrder::MPASEdgeData::print_all()
