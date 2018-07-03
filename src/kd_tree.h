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

//kd tree for quick closest point queries

#ifndef KD_TREE_H_
#define KD_TREE_H_

#include <vector>
#include <float.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <algorithm>

struct kd_point {
	double x,  y;
	double nx, ny;
};

class kd_tree {
private:
	struct kd_node {
		kd_node *left, *right;
		kd_point point;
	};

	double m_bbmin_x = DBL_MAX, m_bbmax_x = -DBL_MAX;
	double m_bbmin_y = DBL_MAX, m_bbmax_y = -DBL_MAX;

	kd_point *m_points = NULL;
	kd_node  *m_root   = NULL;

	kd_node *do_build_tree(kd_point *points, unsigned int start, unsigned int end, unsigned int level);
	void do_find_closest_point(kd_node *cur_node, double qx, double qy, kd_point *closest, double *current_best, unsigned int level);
	void do_wipe_tree(kd_node *node);

public:
	kd_point find_closest_point(double qx, double qy);

	kd_tree(kd_point *points, unsigned int n);
	~kd_tree();
};

#endif /* KD_TREE_H_ */
