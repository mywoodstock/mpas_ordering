/**
 *  Project:
 *
 *  File: mpas_ordering.hpp
 *  Created: Nov 12, 2013
 *  Modified: Fri 13 Dec 2013 01:38:01 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __MPAS_ORDERING_HPP__
#define __MPAS_ORDERING_HPP__

#include <string>
#include <vector>
#include <map>


enum sfc_t {
	NONE,										// no reordering - original
	RANDOM,
	MORTON_SFC,									// morton ordering
	HILBERT_SFC,								// hilbert ordering
	XYZ_SORT									// sort based on x, y, z coords
	// ...
}; // enum sfc_t


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

		typedef MPASElementData mpas_element_t;
		typedef std::vector <mpas_element_t> vec_mpas_element_t;
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
		bool generate_partition_list();
		bool reindex_ordering_index();
		bool reorder_grid(std::string);
		bool save_graph_info(std::string);
		bool save_partition_info(std::string);
		template<typename T> bool reorder_data(T*, T*, int, int, long*);

		bool reorder_elements_morton_sfc();
		bool reorder_elements_hilbert_sfc();
		bool reorder_elements_xyz_sort();
		bool reorder_elements_random();

		// data

		std::string netcdf_grid_filename_;
		std::string graph_info_filename_;
		std::string graph_partition_filename_;

		unsigned int num_cells_;				// number of cells
		unsigned int num_partitions_;			// number of graph partitions
		unsigned int num_edges_;
		vec_mpas_element_t element_list_;		// list of all elements ordered
		map_original_index_t original_index_map_;	// map from current index to current index in list
		map_original_index_t current_index_map_;	// map from original_index_ to current index in list
		partition_map_t partition_list_;		// list of partitions with elements' original index

}; // class MPASElementOrder

#endif // __MPAS_OERDERING_HPP__
