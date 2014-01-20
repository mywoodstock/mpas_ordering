/**
 *  Project:
 *
 *  File: mpas_ordering.hpp
 *  Created: Nov 12, 2013
 *  Modified: Sun 19 Jan 2014 10:17:23 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __MPAS_ORDERING_HPP__
#define __MPAS_ORDERING_HPP__

#include <string>
#include <vector>
#include <map>

#include "mpas_ordering_typedefs.hpp"


class MPASElementOrder {

	public:
		typedef std::vector <unsigned int> vec_uint_t;
		typedef std::map <int, vec_uint_t> partition_map_t;

		class MPASElementCoord {
			public:
				double x_;
				double y_;
				double z_;
				MPASElementCoord(): x_(0.0), y_(0.0), z_(0.0) { }
				MPASElementCoord(double x, double y, double z): x_(x), y_(y), z_(z) { }
				~MPASElementCoord() { }
		}; // class MPASElementCoord

		typedef MPASElementCoord mpas_element_coord_t;
		typedef std::vector <mpas_element_coord_t> vec_mpas_element_coord_t;
		typedef std::pair <mpas_element_coord_t, mpas_element_coord_t> mpas_element_coord_pair_t;

		class MPASElementData {
			public:

			unsigned int partition_num_;		// starts from 0
			unsigned int ordering_index_;		// starts from 0 (ordering within a partition)
												// sorting key is (partition, ordering_index)
			unsigned int original_index_;		// starts from 0

			double x_coord_;
			double y_coord_;
			double z_coord_;
			vec_uint_t neighbor_list_;			// list of neighbors' original indices

			MPASElementData() { }
			~MPASElementData() { }

			bool operator<(const MPASElementData& b) const {
				return (partition_num_ != b.partition_num_) ?
						(partition_num_ < b.partition_num_) : (ordering_index_ < b.ordering_index_);
			} // operator<()

			bool print_all() const;
		}; // class MPASElementData

		class MPASEdgeData {
			public:

			unsigned int partition_num_;	// partition number
			unsigned int ordering_index_;	// current ordering index
			unsigned int original_index_;	// original index in grid file

			unsigned int cell_index_[2];	// cells on this edge (original index as in grid file)

			MPASEdgeData() { }
			~MPASEdgeData() { }

			bool operator<(const MPASEdgeData& b) const {
				return (partition_num_ != b.partition_num_) ?
						(partition_num_ < b.partition_num_) : (ordering_index_ < b.ordering_index_);
			} // operator<()

			bool print_all() const;
		}; // class MPASEdgeData

		// TODO: merge above two classes into one ...

		typedef MPASElementData mpas_element_t;
		typedef std::vector <mpas_element_t> vec_mpas_element_t;
		typedef MPASEdgeData mpas_edge_t;
		typedef std::vector <mpas_edge_t> vec_mpas_edge_t;
		typedef std::map <unsigned int, unsigned int> map_original_index_t;

	public:

		MPASElementOrder(std::string, std::string);					// no partitions
		MPASElementOrder(std::string, std::string, std::string);	// with partitions
		~MPASElementOrder();

		bool reorder_elements(sfc_t);			// perform reordering
		bool save_elements_order(std::string);	// save the current ordering
												// saves the graph and parts data
		bool print_all();

	private:

		bool init();
		bool generate_original_index_map();
		bool generate_edge_original_index_map();
		bool generate_partition_list();
		bool reindex_ordering_index();
		bool reorder_edges();
		bool reorder_grid(std::string);
		bool save_graph_info(std::string);
		bool save_partition_info(std::string);
		template<typename T> bool reorder_data(T*, T*, int, int, long*);
		template<typename T> bool reorder_edge_data(T*, T*, int, int, long*);

		// reordering functions

		bool reorder_elements_random();
		bool reorder_elements_xyz_sort();
		bool reorder_elements_morton_sfc();
		bool reorder_elements_morton_sfc_new();
		bool reorder_elements_hilbert_sfc();
		bool reorder_elements_peano_sfc();

		// helpers

		mpas_element_coord_pair_t minmax(vec_mpas_element_t::const_iterator&,
										vec_mpas_element_t::const_iterator&) const;
		mpas_element_coord_pair_t minmax(vec_mpas_element_coord_t::const_iterator&,
										vec_mpas_element_coord_t::const_iterator&) const;
		bool project_points_to_plane(vec_mpas_element_t::const_iterator&,
										vec_mpas_element_t::const_iterator&,
										mpas_element_coord_pair_t, vec_mpas_element_coord_t&);
		uint64_t decimal_to_ternary(uint64_t);
		uint64_t ternary_to_decimal(uint64_t);
		uint64_t generate_morton_key(unsigned int, unsigned int, unsigned int, unsigned int);
		uint64_t generate_hilbert_key(unsigned int, unsigned int, unsigned int, unsigned int);
		uint64_t generate_peano_key(unsigned int, unsigned int, unsigned int, unsigned int);

		// data

		std::string netcdf_grid_filename_;
		std::string graph_info_filename_;
		std::string graph_partition_filename_;

		unsigned int num_cells_;				// number of cells
		unsigned int num_partitions_;			// number of graph partitions
		unsigned int num_edges_;				// number of edges

		vec_mpas_element_t element_list_;		// list of all elements ordered
		map_original_index_t original_index_map_;	// map from current index to original index in list
		map_original_index_t current_index_map_;	// map from original_index_ to current index in list

		vec_mpas_edge_t edge_list_;				// list of all edge data
		map_original_index_t edge_original_index_map_;	// from current index to original index in list
		map_original_index_t edge_current_index_map_;	// from original_index_ to current index in list

		partition_map_t partition_list_;		// list of partitions with elements' original index

}; // class MPASElementOrder

#endif // __MPAS_OERDERING_HPP__
