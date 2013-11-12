/**
 *  Project:
 *
 *  File: mpas_ordering.hpp
 *  Created: Nov 12, 2013
 *  Modified: Tue 12 Nov 2013 11:00:17 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __MPAS_ORDERING_HPP__
#define __MPAS_ORDERING_HPP__

#include <string>

enum sfc_t {
	NONE,						// no reordering
	MORTON_SFC,					// morton ordering
	HILBERT_SFC,				// hilbert ordering
	XYZ_SORT					// sort based on x, y, z coords
	// ...
}; // enum sfc_t

class MPASElementOrder {

	class MPASElementData {
		unsigned int original_index_;								// starts from 0
		unsigned int ordering_index_;								// starts from 0
																	// (ordering within a partition)
		unsigned int partition_num_;								// starts from 0
		double x_coord_;
		double y_coord_;
		double z_coord_;
		std::vector<unsigned int> neighbor_list_;					// list of original neighbor indices

		MPASElementData() { }
		~MPASElementData() { }
	}; // class MPASElementData

	public:
		MPASElementOrder(std::string, std::string);					// no partitions
		MPASElementOrder(std::string, std::string, std::string);	// with partitions
		~MPASElementOrder();

		bool reorder_elements(sfc_t);								// perform reordering
		bool save_elements_order(std::string);						// save the current ordering
																	// saves the graph and parts data

	private:
		bool init();
		bool reorder_elements_morton_sfc();
		bool reorder_elements_hilbert_sfc();
		bool reorder_elements_xyz_sort();

		// data

		std::string netcdf_grid_filename_;
		std::string graph_info_filename_;
		std::string graph_partition_filename_;

		unsigned int num_cells_;									// number of cells
		unsigned int num_partitions_;								// number of graph partitions
		std::vector<MPASElementData> element_list_;					// list of all elements
		std::map<int, std::vector<unsigned int> > partition_list_;	// partitions

}; // class MPASElementOrder

#endif // __MPAS_OERDERING_HPP__
