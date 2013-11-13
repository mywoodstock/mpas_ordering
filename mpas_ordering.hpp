/**
 *  Project:
 *
 *  File: mpas_ordering.hpp
 *  Created: Nov 12, 2013
 *  Modified: Tue 12 Nov 2013 03:51:34 PM PST
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
	MORTON_SFC,									// morton ordering
	HILBERT_SFC,								// hilbert ordering
	XYZ_SORT									// sort based on x, y, z coords
	// ...
}; // enum sfc_t


class MPASElementOrder {

	public:
		typedef std::vector <unsigned int> vec_uint_t;
		typedef std::map <int, vec_uint_t> partition_map_t;

		class MPASElementData {
			public:

			unsigned int original_index_;			// starts from 0 == index in element_list_
			unsigned int partition_num_;			// starts from 0
			unsigned int ordering_index_;			// starts from 0
												// (ordering within a partition)
			double x_coord_;
			double y_coord_;
			double z_coord_;
			vec_uint_t neighbor_list_;				// list of original neighbor indices

			MPASElementData() { }
			~MPASElementData() { }

			bool operator<(const MPASElementData& b) const {
				return (partition_num_ != b.partition_num_) ?
							(partition_num_ < b.partition_num_) : (ordering_index_ < b.ordering_index_);
			} // operator<()

			bool print_all() const;
		}; // class MPASElementData

		typedef MPASElementData mpas_element_t;
		typedef std::map <unsigned int, mpas_element_t> map_mpas_element_t;
		typedef std::map <unsigned int, unsigned int> map_reverse_index_t;

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
		bool reindex_ordering_index();

		bool reorder_elements_morton_sfc();
		bool reorder_elements_hilbert_sfc();
		bool reorder_elements_xyz_sort();

		// data

		std::string netcdf_grid_filename_;
		std::string graph_info_filename_;
		std::string graph_partition_filename_;

		unsigned int num_cells_;				// number of cells
		unsigned int num_partitions_;			// number of graph partitions
		map_mpas_element_t element_list_;		// list of all elements ordered with "map_index"
		map_reverse_index_t original_index_map_;	// map from "map_index" to original_index_
		partition_map_t partition_list_;		// partitions with sorted element indices

}; // class MPASElementOrder

#endif // __MPAS_OERDERING_HPP__
