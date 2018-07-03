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

#include "element.h"

element::element(unsigned int ll_idx, unsigned int lr_idx, unsigned int tr_idx, unsigned int tl_idx) :
		ll_idx(ll_idx), lr_idx(lr_idx), tr_idx(tr_idx), tl_idx(tl_idx) {}

std::vector<element> element_map_rectangle(unsigned int ex, unsigned int ey) {
	std::vector<element> map;

	unsigned int nx = ex+1;
	const unsigned int pattern[4] = {0, 1, nx+1, nx};
	unsigned int inc = 0;

	for (unsigned int i = 0; i < ey; i++) {
		for (unsigned int j = 0; j < ex; j++) {
			map.push_back(element(pattern[0]+inc,pattern[1]+inc,pattern[2]+inc,pattern[3]+inc));
			inc++;
		}
		inc = (i+1)*nx;
	}

	return map;
}

