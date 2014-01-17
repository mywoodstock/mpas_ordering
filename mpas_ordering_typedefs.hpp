/**
 *  Project:
 *
 *  File: mpas_ordering_typedefs.hpp
 *  Created: Jan 16, 2014
 *  Modified: Thu 16 Jan 2014 11:22:04 AM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#ifndef __MPAS_ORDERING_TYPEDEFS_HPP__
#define __MPAS_ORDERING_TYPEDEFS_HPP__


enum sfc_t {
	NONE,										// no reordering - original
	RANDOM,										// completely randomized
	MORTON_SFC,									// morton ordering
	HILBERT_SFC,								// hilbert ordering
	XYZ_SORT									// sort based on x, y, z coords
	// ...
}; // enum sfc_t


#endif // __MPAS_ORDERING_TYPEDEFS_HPP__
