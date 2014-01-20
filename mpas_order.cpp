/**
 *  Project:
 *
 *  File: mpas_order.cpp
 *  Created: Nov 12, 2013
 *  Modified: Mon 20 Jan 2014 09:06:45 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <string>

#include "mpas_ordering.hpp"

int main(int narg, char** args) {
	char *part_filename;
	std::string ordering;

	if(narg == 6) {
		part_filename = args[4];
		ordering = args[5];
	} else if(narg == 5) {
		//part_filename = (char*) "";
		std::cout << "usage: mpas_order <grid_file> <graph_file> <output_prefix> "
					<< "<partition_file> <ordering>" << std::endl;
		return 1;
	} else {
		std::cout << "usage: mpas_order <grid_file> <graph_file> <output_prefix> "
					<< "<partition_file> <ordering>" << std::endl;
		return 1;
	} // if-else

	MPASElementOrder my_ordering(args[1], args[2], part_filename);
	sfc_t type;
	if(ordering.compare("morton") == 0)			type = MORTON_SFC;
	else if(ordering.compare("random") == 0)	type = RANDOM;
	else if(ordering.compare("hilbert") == 0)	type = HILBERT_SFC;
	else if(ordering.compare("peano") == 0)		type = PEANO_SFC;
	else if(ordering.compare("xyz") == 0)		type = XYZ_SORT;
	else {
		std::cerr << "error: invalid ordering given. choices are "
					<< "'xyz', 'random', 'hilbert', 'morton', 'peano'" << std::endl;
		return 1;
	} // if-else
	my_ordering.reorder_elements(type);
	my_ordering.save_elements_order(args[3]);

	return 0;
} // main()
