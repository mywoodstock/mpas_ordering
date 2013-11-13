/**
 *  Project:
 *
 *  File: mpas_order.cpp
 *  Created: Nov 12, 2013
 *  Modified: Wed 13 Nov 2013 11:45:32 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>

#include "mpas_ordering.hpp"

int main(int narg, char** args) {
	char* part_filename;

	if(narg == 4) {
		part_filename = args[3];
	} else if(narg == 3) {
		part_filename = "";
	} else {
		std::cout << "usage: mpas_order <grid_file> <graph_file> [<partition_file>]" << std::endl;
		return 1;
	} // if-else

	MPASElementOrder my_ordering(args[1], args[2], part_filename);
	my_ordering.print_all();
	my_ordering.save_elements_order("hehehehe.info");

	return 0;
} // main()
