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

#ifndef UTILS_H_
#define UTILS_H_

template<typename T>
T** mesh_ring(unsigned int nbox, unsigned int &n, std::vector<element> &elements_in, double rin, double rout) {

	double dx = 2.*rout/(nbox-1);

	std::vector<std::tuple<particle, std::vector<unsigned int>>> points(nbox*nbox);

	// create points
	unsigned int part_iter = 0;
	for(unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			std::get<0>(points[part_iter]).idx = part_iter;
			std::get<0>(points[part_iter]).x = i*dx - rout;
			std::get<0>(points[part_iter]).y = j*dx - rout;
			part_iter++;
		}
	}

	// element map for those points
	std::vector<element> elements = element_map_rectangle(nbox-1, nbox-1);

	// compact element map to contain only elements which are inside the ring
	for (auto it = elements.begin(); it != elements.end(); ++it) {
		particle ll = std::get<0>(points[it->ll_idx]);
		particle lr = std::get<0>(points[it->lr_idx]);
		particle tr = std::get<0>(points[it->tr_idx]);
		particle tl = std::get<0>(points[it->tl_idx]);

		double r_ll = ll.x*ll.x + ll.y*ll.y;
		double r_lr = lr.x*lr.x + lr.y*lr.y;
		double r_tr = tr.x*tr.x + tr.y*tr.y;
		double r_tl = tl.x*tl.x + tl.y*tl.y;

		bool in_ll = r_ll < rout*rout && r_ll > rin*rin;
		bool in_lr = r_lr < rout*rout && r_lr > rin*rin;
		bool in_tr = r_tr < rout*rout && r_tr > rin*rin;
		bool in_tl = r_tl < rout*rout && r_tl > rin*rin;

		if (in_ll && in_lr && in_tr && in_tl) {
			elements_in.push_back(*it);
		}
	}

	// install back pointers on the particles
	for (unsigned int i = 0; i < elements_in.size(); i++) {
		unsigned int ll_idx = elements_in[i].ll_idx;
		unsigned int lr_idx = elements_in[i].lr_idx;
		unsigned int tr_idx = elements_in[i].tr_idx;
		unsigned int tl_idx = elements_in[i].tl_idx;

		std::get<1>(points[ll_idx]).push_back(i);
		std::get<1>(points[lr_idx]).push_back(i);
		std::get<1>(points[tr_idx]).push_back(i);
		std::get<1>(points[tl_idx]).push_back(i);
	}

	// compact points with awesume c++ 11 magic
	std::vector<std::tuple<particle, std::vector<unsigned int>>> points_in(nbox*nbox);
	auto it = std::copy_if(points.begin(), points.end(), points_in.begin(),
			[](std::tuple<particle, std::vector<unsigned int>> t) { return std::get<1>(t).size() > 0;});
	points_in.resize(std::distance(points_in.begin(),it));

	n = points_in.size();

	// follow the back pointers and fix them to reflect the new positions
	for (unsigned int i = 0; i < points_in.size(); i++) {
		unsigned int new_idx = i;
		unsigned int old_idx = std::get<0>(points_in[i]).idx;

		std::vector<unsigned int> back_pointers = std::get<1>(points_in[i]);
		for (auto it = back_pointers.begin(); it != back_pointers.end(); ++it) {
			if (elements_in[*it].ll_idx == old_idx) elements_in[*it].ll_idx = new_idx;
			if (elements_in[*it].lr_idx == old_idx) elements_in[*it].lr_idx = new_idx;
			if (elements_in[*it].tr_idx == old_idx) elements_in[*it].tr_idx = new_idx;
			if (elements_in[*it].tl_idx == old_idx) elements_in[*it].tl_idx = new_idx;
		}
	}

	T **particles = new T*[n];
	for (unsigned int i = 0; i < n; i++) {
		particles[i] = new T(i);
		particles[i]->x = std::get<0>(points_in[i]).x;
		particles[i]->y = std::get<0>(points_in[i]).y;
		particles[i]->X = std::get<0>(points_in[i]).x;
		particles[i]->Y = std::get<0>(points_in[i]).y;
	}

	return particles;
}

#endif /* UTILS_H_ */
