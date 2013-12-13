/**
 *  Project:
 *
 *  File: mpas_ordering_random.cpp
 *  Created: Dec 13, 2013
 *  Modified: Fri 13 Dec 2013 01:33:53 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>

#include "mpas_ordering.hpp"


bool MPASElementOrder::reorder_elements_random() {
	// randomize
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::shuffle(element_list_.begin(), element_list_.end(), std::default_random_engine(seed));
	// reindex and sort according to partitions
	reindex_ordering_index();
	std::sort(element_list_.begin(), element_list_.end());
	return true;
} // MPASElementOrder::reorder_elements_random()
