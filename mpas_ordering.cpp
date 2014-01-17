/**
 *  Project:
 *
 *  File: mpas_ordering.cpp
 *  Created: Nov 12, 2013
 *  Modified: Thu 16 Jan 2014 05:51:59 PM PST
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


//bool comp_umap_elements(const std::pair<int, MPASElementOrder::mpas_element_t>& a,
//					const std::pair<int, MPASElementOrder::mpas_element_t>& b) {
//	return (a.second < b.second);
//} // comp_umap_elements()


bool MPASElementOrder::reorder_elements(sfc_t sfc) {
	switch(sfc) {
		case RANDOM:
			return reorder_elements_random();

		case MORTON_SFC:
			return reorder_elements_morton_sfc_new();

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
	// create the mappers
	generate_original_index_map();

#ifdef DEBUG
	for(int i = 0; i < element_list_.size(); ++ i)
		std::cout << original_index_map_[i] + 1 << " -> " << i + 1 << std::endl;
#endif

	// reorder stuff in the grid and write a new grid file
	reorder_grid(filename_prefix);
	save_graph_info(filename_prefix);
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

	num_edges_ = 0;
	for(int i = 0; i < num_cells_; ++ i) {
		MPASElementData cell;
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
		num_edges_ += cell.neighbor_list_.size();
		element_list_.push_back(cell);
	} // for
	num_edges_ = num_edges_ >> 1;

	delete[] z_cells;
	delete[] y_cells;
	delete[] x_cells;

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


bool MPASElementOrder::reorder_grid(std::string prefix) {
	std::stringstream grid_name;
	if(num_partitions_ > 1)
		grid_name << prefix << "." << num_partitions_ << ".nc";
	else
		grid_name << prefix << ".nc";
	#ifdef _64BITOPFFSET
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

	// read vars and reorder data
	NcDim *dims[5];
	long dim_sizes[5];
	for(int i = 0; i < ncid.num_vars(); ++ i) {
		NcVar *var_id = ncid.get_var(i);
		int num_dims = (*var_id).num_dims();
		int size = 1;
		std::map <NcToken, int> var_dims;
		bool to_reorder = false;
		int cell_dim_num = -1;
		for(int j = 0; j < num_dims; ++ j) {
			dims[j] = (*var_id).get_dim(j);
			dim_sizes[j] = (*dims[j]).size();
			size *= (*dims[j]).size();
			var_dims[(*dims[j]).name()] = (*dims[j]).size();
			if((*dims[j]).name() == cell_dim_token && (*var_id).name() != cell_id_token) {
				to_reorder = true;
				cell_dim_num = j;
			} // if
		} // for
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
				if(to_reorder) {
					int_buff_out = new (std::nothrow) int[size];
					reorder_data<int>(int_buff_in, int_buff_out, cell_dim_num, num_dims, dim_sizes);
				} else {
					int_buff_out = int_buff_in;
				} // if-else
				// update the cell numbers in cellsOnEdge, cellsOnCell, cellsOnVertex
				if((*var_id).name() == cells_on_edge_token ||
						(*var_id).name() == cells_on_cell_token ||
						(*var_id).name() == cells_on_vertex_token) {
					for(unsigned int i = 0; i < size; ++ i)
						if(int_buff_out[i] > 0)
							int_buff_out[i] = current_index_map_[int_buff_out[i] - 1] + 1;
				} // for
				(*out_var_id).put(int_buff_out, dim_sizes);
				if(to_reorder) delete[] int_buff_out;
				delete[] int_buff_in;
				break;

			case ncDouble:
				dbl_buff_in = new (std::nothrow) double[size];
				(*var_id).get(dbl_buff_in, dim_sizes);
				if(to_reorder) {
					dbl_buff_out = new (std::nothrow) double[size];
					reorder_data<double>(dbl_buff_in, dbl_buff_out, cell_dim_num, num_dims, dim_sizes);
					(*out_var_id).put(dbl_buff_out, dim_sizes);
					delete[] dbl_buff_out;
				} else {
					(*out_var_id).put(dbl_buff_in, dim_sizes);
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
