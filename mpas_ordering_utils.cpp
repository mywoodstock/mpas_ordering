/**
 *  Project:
 *
 *  File: mpas_ordering_utils.cpp
 *  Created: Jan 16, 2014
 *  Modified: Sat 18 Jan 2014 04:46:59 PM PST
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <limits>
#include <algorithm>
#include <cmath>

#include "mpas_ordering.hpp"


typedef MPASElementOrder::mpas_element_coord_pair_t mpas_element_coord_pair_t;

// for a given range in element_list_, find min and max coordinates
mpas_element_coord_pair_t
MPASElementOrder::minmax(vec_mpas_element_t::const_iterator& begin,
						vec_mpas_element_t::const_iterator& end) const {
	double d_max = std::numeric_limits<double>::max();
	double d_min = std::numeric_limits<double>::min();
	mpas_element_coord_t min_point(d_max, d_max, d_max), max_point(d_min, d_min, d_min);
	for(vec_mpas_element_t::const_iterator i = begin; i != end; ++ i) {
		min_point.x_ = (min_point.x_ > (*i).x_coord_) ? (*i).x_coord_ : min_point.x_;
		max_point.x_ = (max_point.x_ < (*i).x_coord_) ? (*i).x_coord_ : max_point.x_;
		min_point.y_ = (min_point.y_ > (*i).y_coord_) ? (*i).y_coord_ : min_point.y_;
		max_point.y_ = (max_point.y_ < (*i).y_coord_) ? (*i).y_coord_ : max_point.y_;
		min_point.z_ = (min_point.z_ > (*i).z_coord_) ? (*i).z_coord_ : min_point.z_;
		max_point.z_ = (max_point.z_ < (*i).z_coord_) ? (*i).z_coord_ : max_point.z_;
	} // for
	return mpas_element_coord_pair_t(min_point, max_point);
} // MPASElementOrder::minmax()


// for a given range of coords, find min and max
mpas_element_coord_pair_t
MPASElementOrder::minmax(vec_mpas_element_coord_t::const_iterator& begin,
						vec_mpas_element_coord_t::const_iterator& end) const {
	vec_mpas_element_coord_t::const_iterator i = begin;
	mpas_element_coord_t min_point((*i).x_, (*i).y_, (*i).z_), max_point((*i).x_, (*i).y_, (*i).z_);
	++ i;
	for(; i != end; ++ i) {
		min_point.x_ = (min_point.x_ > (*i).x_) ? (*i).x_ : min_point.x_;
		max_point.x_ = (max_point.x_ < (*i).x_) ? (*i).x_ : max_point.x_;
		min_point.y_ = (min_point.y_ > (*i).y_) ? (*i).y_ : min_point.y_;
		max_point.y_ = (max_point.y_ < (*i).y_) ? (*i).y_ : max_point.y_;
		min_point.z_ = (min_point.z_ > (*i).z_) ? (*i).z_ : min_point.z_;
		max_point.z_ = (max_point.z_ < (*i).z_) ? (*i).z_ : max_point.z_;
	} // for
	return mpas_element_coord_pair_t(min_point, max_point);
} // MPASElementOrder::minmax()


bool pair_compare(std::pair <unsigned int, unsigned int> a, std::pair <unsigned int, unsigned int> b) {
	return a.first < b.first;
} // pair_compare


bool MPASElementOrder::project_points_to_plane(vec_mpas_element_t::const_iterator& begin,
							vec_mpas_element_t::const_iterator& end,
							mpas_element_coord_pair_t minmax_coords,
							vec_mpas_element_coord_t& projected_coords) {
	projected_coords.clear();
	mpas_element_coord_t center_point, projection;
	center_point.x_ = (minmax_coords.first.x_ + minmax_coords.second.x_) / 2;
	center_point.y_ = (minmax_coords.first.y_ + minmax_coords.second.y_) / 2;
	center_point.z_ = (minmax_coords.first.z_ + minmax_coords.second.z_) / 2;
	//std::cout << "min: (" << minmax_coords.first.x_ << ", "
	//					<< minmax_coords.first.y_ << ", "
	//					<< minmax_coords.first.z_ << ") max: ("
	//					<< minmax_coords.second.x_ << ", "
	//					<< minmax_coords.second.y_ << ", "
	//					<< minmax_coords.second.z_ << ") center: ("
	//					<< center_point.x_ << ", "
	//					<< center_point.y_ << ", "
	//					<< center_point.z_ << ")" << std::endl;
	double norm = center_point.x_ * center_point.x_ +
					center_point.y_ * center_point.y_ +
					center_point.z_ * center_point.z_;
	double a = center_point.x_ / norm;
	double b = center_point.y_ / norm;
	double c = center_point.z_ / norm;
	double abc = a * a + b * b + c * c;
	double k = c / sqrt(abc);
	double s = sqrt(1.0 - k * k);
	double l = 1.0 - k;
	double X = b / sqrt(a * a + b * b);
	double Y = - a / sqrt(a * a + b * b);
	double Z = 0.0;
	// rotation matrix
	double r00 = X * X * l + k,     r01 = X * Y * l - Z * s, r02 = X * Z * l + Y * s;
	double r10 = Y * X * l + Z * s, r11 = Y * Y * l + k,     r12 = Y * Z * l - X * s;
	double r20 = Z * X * l - Y * s, r21 = Z * Y * l + X * s, r22 = Z * Z * l + k;
	for(vec_mpas_element_t::const_iterator i = begin; i != end; ++ i) {
		double xp = (*i).x_coord_, yp = (*i).y_coord_, zp = (*i).z_coord_;
		double t = (1.0 - (a * xp + b * yp + c * zp)) / abc;
		double x1 = xp + a * t;
		double y1 = yp + b * t;
		double z1 = zp + c * t;
		projection.x_ = r00 * x1 + r01 * y1 + r02 * z1;
		projection.y_ = r10 * x1 + r11 * y1 + r12 * z1;
		projection.z_ = 0.0; //r20 * x1 + r21 * y1 + r22 * z1;
		projected_coords.push_back(projection);
		//std::cout << "*@@ [ " << xp << " , " << yp << " , " << zp << " ] -> [ "
		//			<< projection.x_ << " , " << projection.y_ << " , " << projection.z_ << " ]"
		//			<< std::endl;
	} // for

//	// make them all positive
//	vec_mpas_element_coord_t::const_iterator cb = projected_coords.begin();
//	vec_mpas_element_coord_t::const_iterator ce = projected_coords.end();
//	mpas_element_coord_pair_t new_minmax = minmax(cb, ce);
//	for(vec_mpas_element_coord_t::iterator i = projected_coords.begin(); i != projected_coords.end(); ++ i) {
//		if(new_minmax.first.x_ < 0.0) (*i).x_ = (*i).x_ - new_minmax.first.x_;
//		if(new_minmax.first.y_ < 0.0) (*i).y_ = (*i).y_ - new_minmax.first.y_;
//		//if(new_minmax.first.z_ < 0.0) (*i).z_ = (*i).z_ - new_minmax.first.z_;
//	} // for

	// sort on x and re-assign coords
	// sort on y and re-assign coords
	std::vector <std::pair <unsigned int, unsigned int> > x_mapper, y_mapper;
	unsigned int idx = 0;
	for(vec_mpas_element_coord_t::iterator i = projected_coords.begin();
		   	i != projected_coords.end(); ++ i, ++ idx) {
		x_mapper.push_back(std::pair <unsigned int, unsigned int> ((*i).x_, idx));
		y_mapper.push_back(std::pair <unsigned int, unsigned int> ((*i).y_, idx));
	} // for
	std::sort(x_mapper.begin(), x_mapper.end(), pair_compare);
	std::sort(y_mapper.begin(), y_mapper.end(), pair_compare);
	for(idx = 0; idx < projected_coords.size(); ++ idx) {
		projected_coords[x_mapper[idx].second].x_ = idx;
		projected_coords[y_mapper[idx].second].y_ = idx;
	} // for

	return true;
} // MPASElementOrder::project_points_to_plane()


uint64_t MPASElementOrder::decimal_to_ternary(uint64_t n) {
	uint64_t t = 0, ten = 1;
	for(int i = 0; i < 64 && n != 0; ++ i) {
		t += (n % 3) * ten;
		n = n / 3;
		ten *= 10;
	} // for
	return t;
} // MPASElementOrder::decimal_to_ternary()


uint64_t MPASElementOrder::ternary_to_decimal(uint64_t n) {
	uint64_t d = 0, three = 1;
	for(int i = 0; i < 64; ++ i) {
		d += (n % 10) * three;
		three *= 3;
		n /= 10;
	} // for
	return d;
} // MPASElementOrder::decimal_to_ternary()
