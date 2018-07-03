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

#include "kd_tree.h"

kd_tree::kd_tree(kd_point *points, unsigned int n) {
	for (unsigned int i = 0; i < n; i++) {
		m_bbmin_x = fmin(m_bbmin_x, points[i].x);
		m_bbmin_y = fmin(m_bbmin_y, points[i].y);
		m_bbmax_x = fmax(m_bbmax_x, points[i].x);
		m_bbmax_y = fmax(m_bbmax_y, points[i].y);
	}

	m_points = new kd_point[n];
	for (unsigned int i = 0; i < n; i++) {
		m_points[i] = points[i];
	}

	m_root = do_build_tree(m_points, 0, n, 0);
}

kd_tree::kd_node *kd_tree::do_build_tree(kd_point *points, unsigned int start, unsigned int end, unsigned int level) {
	kd_node *node = new kd_node();

	std::sort(points + start, points + end, (level % 2) ? [](const kd_point a, const kd_point b) {return a.x < b.x;} :
			[](const kd_point a, const kd_point b) {return a.y < b.y;});

	 int median = (end+start)/2;

	node->point = points[median];

	node->left  = NULL;
	node->right = NULL;

	int left_lo = start;
	int left_hi = median;

	int right_lo = median+1;
	int right_hi = end;

	node->left  = (left_hi  - left_lo  > 0) ? do_build_tree(points, left_lo,  left_hi,  level+1) : NULL;
	node->right = (right_hi - right_lo > 0) ? do_build_tree(points, right_lo, right_hi, level+1) : NULL;

	return node;
}

void kd_tree::do_find_closest_point(kd_node *cur_node, double qx, double qy, kd_point *closest, double *current_best, unsigned int level) {
	if (!cur_node) return;

	bool left = (level % 2 && qx < cur_node->point.x) || (!(level % 2) && qy < cur_node->point.y);

	if (left) {
		do_find_closest_point(cur_node->left,  qx, qy, closest, current_best, level+1);
	} else {
		do_find_closest_point(cur_node->right, qx, qy, closest, current_best, level+1);
	}

	//unwind
	double dx = qx - cur_node->point.x;
	double dy = qy - cur_node->point.y;

	double dx2 = dx*dx;
	double dy2 = dy*dy;

	double dist2 = dx2 + dy2;

	if (dist2 < *current_best) {
		*current_best = dist2;
		*closest = cur_node->point;
	}

	//if we are too close to a a cutting plane we need to look at the other side since there might be a closer point
	if (*current_best > dx2 || *current_best > dy2) {
		if (left) {
			do_find_closest_point(cur_node->right, qx, qy, closest, current_best, level+1);
		} else {
			do_find_closest_point(cur_node->left,  qx, qy, closest, current_best, level+1);
		}
	}
}

kd_point kd_tree::find_closest_point(double qx, double qy) {
	kd_point closest;
	double best = DBL_MAX;
	do_find_closest_point(m_root, qx, qy, &closest, &best, 0);
	return closest;
}

void kd_tree::do_wipe_tree(kd_node *node) {
	if (!node) return;

	do_wipe_tree(node->left);
	do_wipe_tree(node->right);

	if (node->left)  delete node->left;
	if (node->right) delete node->right;
}

kd_tree::~kd_tree() {
	do_wipe_tree(m_root);
	delete[] m_points;
}
