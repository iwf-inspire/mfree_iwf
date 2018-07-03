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

#ifndef RIEMANN_INTEGRATION_H_
#define RIEMANN_INTEGRATION_H_

template<typename T>
void riemann_update_integration_weights_gp(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		b.get_cur_quad_points()[i]->quad_weight = b.get_cur_quad_points()[i]->m/b.get_cur_quad_points()[i]->rho;
	}
}

template<typename T>
void riemann_update_integration_weights(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		b.get_cur_particles()[i]->quad_weight = b.get_cur_particles()[i]->m/b.get_cur_particles()[i]->rho;
	}
}

#endif /* RIEMANN_INTEGRATION_H_ */
