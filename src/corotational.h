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

//code to compute the shape matching matrix described in
//Becker, Markus, Markus Ihmsen, and Matthias Teschner. "Corotated SPH for Deformable Solids." NPH. 2009. | Equation (12)

#ifndef COROTATIONAL_H_
#define COROTATIONAL_H_

#include <math.h>
#include <stdio.h>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
void corotational_anisotropy_tensor(body<T> &b) {
	static_assert(std::is_base_of<particle_corotational, T>::value, "T ist not derived from particle_corotational");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Xi = pi->X;
		double Yi = pi->Y;

		double xi = pi->x;
		double yi = pi->y;

		double Axx = 0.;
		double Axy = 0.;
		double Ayx = 0.;
		double Ayy = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double Xj = pj->X;
			double Yj = pj->Y;

			double xj = pj->x;
			double yj = pj->y;

			double mj = pj->m;

			double xji = xj-xi;
			double yji = yj-yi;

			double Xji = Xj-Xi;
			double Yji = Yj-Yi;

			Axx += xji*Xji*mj*w.w;
			Axy += xji*Yji*mj*w.w;
			Ayx += yji*Xji*mj*w.w;
			Ayy += yji*Yji*mj*w.w;
		}

		double x = Axx + Ayy;
		double y = Ayx - Axy;

		double theta = atan2(y,x);

		pi->Rxx =  cos(theta);
		pi->Rxy = -sin(theta);
		pi->Ryx =  sin(theta);
		pi->Ryy =  cos(theta);
	}
}

#endif /* COROTATIONAL_H_ */
