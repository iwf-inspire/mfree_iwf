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

//interpolate velocity from particles to quadrature points

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
void interpolation_interpolate_velocity(body<T> &b) {
	static_assert(std::is_base_of<particle_ul_weak, T>::value, "T ist not derived from particle_ul_weak");

	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double vxi = 0.;
		double vyi = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj = pj->vx;
			double vyj = pj->vy;

			double quad_weight = pj->quad_weight;

			vxi += vxj*w.w*quad_weight;
			vyi += vyj*w.w*quad_weight;
		}

		pi->vx = vxi;
		pi->vy = vyi;
	}
}

#endif /* INTERPOLATION_H_ */
