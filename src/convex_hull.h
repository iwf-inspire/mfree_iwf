//Copyright ETH Zurich, IWF

//This file is part of mfree_iwf.

//mfree_iwf is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//mfree_iwf/ul_cut is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with mfree_iwf.  If not, see <http://www.gnu.org/licenses/>.

//fast convex hull computation using the monotone chain algorithm
//	adapted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
//	license unclear

#ifndef CONVEX_HULL_H_
#define CONVEX_HULL_H_

#include <vector>
#include <tuple>
#include <algorithm>

static double ccw(std::tuple<double, double, unsigned int>  p1,
		std::tuple<double, double, unsigned int>  p2,
		std::tuple<double, double, unsigned int>  p3) {

	double p1_x = std::get<0>(p1);
	double p1_y = std::get<1>(p1);

	double p2_x = std::get<0>(p2);
	double p2_y = std::get<1>(p2);

	double p3_x = std::get<0>(p3);
	double p3_y = std::get<1>(p3);

	return (p2_x - p1_x)*(p3_y - p1_y) - (p2_y - p1_y)*(p3_x - p1_x);
}

static bool leftlower(std::tuple<double, double, unsigned int> p1, std::tuple<double, double, unsigned int> p2)
{
	double p1_x = std::get<0>(p1);
	double p1_y = std::get<1>(p1);

	double p2_x = std::get<0>(p2);
	double p2_y = std::get<1>(p2);

	if (p1_x < p2_x) return true;
	if (p1_x > p2_x) return false;

	if (p1_y < p2_y) return true;
	if (p1_y > p2_y) return false;

	return 0;
}

template<typename T>
std::vector<unsigned int> convex_hull(T **points, unsigned int num_points) {
	std::vector<std::tuple<double,double,unsigned int>> m_points(num_points);

	for (int i = 0; i < (int) num_points; i++) {
		m_points[i] = std::tuple<double, double, unsigned int>(points[i]->x, points[i]->y, i);
	}

	std::sort(m_points.begin(), m_points.end(), &leftlower);

	std::vector<std::tuple<double, double, unsigned int>> hull(num_points);
	int k = 0;

	//upper hull
	for (int i = 0; i < (int) num_points; ++i) {
		while (k >= 2 && ccw(hull[k-2], hull[k-1], m_points[i]) < 0)
			--k;
		hull[k++] = m_points[i];
	}

	//lower hull
	for (int i = (int) num_points-2, t = k+1; i >= 0; --i) {
		while (k >= t && ccw(hull[k-2], hull[k-1], m_points[i]) < 0)
			--k;
		hull[k++] = m_points[i];
	}

	std::vector<unsigned int> hull_idx;
	for (int i = 0; i < k-1; i++) {
		hull_idx.push_back(std::get<2>(hull[i]));
	}

	return hull_idx;
}

#endif /* CONVEX_HULL_H_ */
