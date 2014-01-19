/**
 *  Project:
 *
 *  File: mpas_ordering_typedefs.hpp
 *  Created: Jan 16, 2014
 *  Modified: Sat 18 Jan 2014 05:23:12 PM PST
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
	PEANO_SFC,									// peano curve ordering
	HILBERT_MOORE_SFC,
	PEANO_MEANDER_SFC,
	SIERPINSKI_SFC,
	LEBESGUE_SFC,
	H_INDEX,
	BETA_OMEGA_SFC,
	GOSPER_FLOWSNAKE_SFC,
	XYZ_SORT									// sort based on x, y, z coords
	// ...
}; // enum sfc_t


#endif // __MPAS_ORDERING_TYPEDEFS_HPP__
