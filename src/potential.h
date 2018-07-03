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

//methods to implement
// Müller, Matthias, et al. "Point based animation of elastic, plastic and melting objects." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium on Computer animation. Eurographics Association, 2004.

#ifndef POTENTIAL_H_
#define POTENTIAL_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
void potential_comp_elastic_potential(body<T> &b) {
	static_assert(std::is_base_of<particle_potential, T>::value, "T ist not derived from particle_potential");


	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		glm::dmat2x2 Ji(pi->Hxx+1., pi->Hyx, pi->Hxy, pi->Hyy+1.);
		glm::dmat2x2 Si(pi->Sxx, pi->Sxy, pi->Sxy, pi->Syy);

		double voli = pi->quad_weight;

		glm::dmat2x2 Fei = -2*voli*Ji*Si;

		double vx_t = 0.;
		double vy_t = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];

			kernel_result wij = pi->w[j];
			kernel_result wji = pi->wji[j];

			T *pj = b.get_cur_particles()[jdx];

			glm::dvec2 grad_wij(wij.w_x, wij.w_y);
			glm::dvec2 fij = Fei*grad_wij;

			glm::dmat2x2 Jj(pj->Hxx+1., pj->Hyx, pj->Hxy, pj->Hyy+1.);
			glm::dmat2x2 Sj(pj->Sxx, pj->Sxy, pj->Sxy, pj->Syy);

			double volj = pj->quad_weight;

			glm::dmat2x2 Fej = -2*volj*Jj*Sj;

			glm::dvec2 grad_wji(wji.w_x, wji.w_y);
			glm::dvec2 fji = Fej*grad_wji;

			vx_t += 0.5*(fji[0] - fij[0])/pi->m;
			vy_t += 0.5*(fji[1] - fij[1])/pi->m;
		}

		pi->vx_t = vx_t;
		pi->vy_t = vy_t;
	}

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];
		pi->vx_t += pi->fcx/pi->m;
		pi->vy_t += pi->fcy/pi->m;
	}
}

#endif /* POTENTIAL_H_ */
