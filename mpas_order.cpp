/**
 *  Project:
 *
 *  File: mpas_order.cpp
 *  Created: Nov 12, 2013
 *  Modified: Fri 13 Dec 2013 10:23:10 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>

#include "mpas_ordering.hpp"

int main(int narg, char** args) {
	char* part_filename;

	if(narg == 5) {
		part_filename = args[4];
	} else if(narg == 4) {
		part_filename = "";
	} else {
		std::cout << "usage: mpas_order <grid_file> <graph_file> <output_prefix> [<partition_file>]" << std::endl;
		return 1;
	} // if-else

	MPASElementOrder my_ordering(args[1], args[2], part_filename);
	my_ordering.reorder_elements(XYZ_SORT);
	//my_ordering.print_all();
	my_ordering.save_elements_order(args[3]);

	return 0;
} // main()
