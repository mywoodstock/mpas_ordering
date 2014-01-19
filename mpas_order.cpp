/**
 *  Project:
 *
 *  File: mpas_order.cpp
 *  Created: Nov 12, 2013
 *  Modified: Sat 18 Jan 2014 04:51:46 PM PST
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
		//part_filename = (char*) "";
		std::cout << "usage: mpas_order <grid_file> <graph_file> <output_prefix> <partition_file>" << std::endl;
		return 1;
	} else {
		std::cout << "usage: mpas_order <grid_file> <graph_file> <output_prefix> <partition_file>" << std::endl;
		return 1;
	} // if-else

	MPASElementOrder my_ordering(args[1], args[2], part_filename);
	//my_ordering.reorder_elements(MORTON_SFC);
	//my_ordering.reorder_elements(RANDOM);
	//my_ordering.reorder_elements(HILBERT_SFC);
	my_ordering.reorder_elements(PEANO_SFC);
	//my_ordering.print_all();
	my_ordering.save_elements_order(args[3]);

	return 0;
} // main()
