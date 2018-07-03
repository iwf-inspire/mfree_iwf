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

//implements the neighbor ellipse idea presented in the paper
//	section 2.2

#ifndef NEIGHBOR_ELLIPSE_H_
#define NEIGHBOR_ELLIPSE_H_

#include <assert.h>
#include <glm/glm.hpp>

#include "particle.h"

static void eigen(glm::dmat2x2 A, glm::dmat2x2 &V, glm::dmat2x2 &D) {

	//diagonal and "almost" diagonal case
	if (fabs(A[1][0]) < 1e-9 && fabs(A[0][1]) < 1e-9) {
		V = glm::dmat2x2(1.); //identity;
		D = glm::dmat2x2(A[0][0], 0., 0., A[1][1]);
		return;
	}

	double T = A[0][0]+A[1][1];
	double det = glm::determinant(A);

	double L1 = T/2 + sqrt(T*T/4-det);
	double L2 = T/2 - sqrt(T*T/4-det);

	double lambda1 = fmin(L1, L2);
	double lambda2 = fmax(L1, L2);

	D[0][0] = lambda1;
	D[1][1] = lambda2;

	if (fabs(A[1][0]) > 1e-12) {
		glm::dvec2 v1 = glm::normalize(glm::dvec2(lambda1 - A[1][1], A[0][1]));
		glm::dvec2 v2 = glm::normalize(glm::dvec2(lambda2 - A[1][1], A[0][1]));
		V = glm::dmat2x2(v1,v2);
	} else { //if (fabs(A[0][1]) > 1e-12)
		glm::dvec2 v1 = glm::normalize(glm::dvec2(A[1][0], lambda1 - A[0][0]));
		glm::dvec2 v2 = glm::normalize(glm::dvec2(A[1][0], lambda2 - A[0][0]));
		V = glm::dmat2x2(v1,v2);
	}
}

template<typename T>
void nbh_ellipse(T **p, unsigned int np, T **q) {
	for (unsigned int i = 0; i < np; i++) {
		T *qp = p[i];

		double px = qp->x;
		double py = qp->y;

		double Hxx = 0.; double Hxy = 0.;
		double Hyx = 0.; double Hyy = 0.;

		double h = 0.;

		for (unsigned int j = 0; j < qp->num_cand_nbh; j++) {
			unsigned int jdx = qp->cand_nbh[j];

			double qx = q[jdx]->x;
			double qy = q[jdx]->y;

			double dx = px - qx;
			double dy = py - qy;

			double r2 = dx*dx + dy*dy;

			h = fmax(h, r2);

			Hxx += dx*dx; Hxy += dx*dy;
			Hyx += dy*dx; Hyy += dy*dy;
		}

		h = sqrt(h);
		h /=2;
		qp->h = h;

		glm::dmat2x2 V,D;
		eigen(glm::dmat2x2(Hxx/qp->num_cand_nbh, Hyx/qp->num_cand_nbh, Hxy/qp->num_cand_nbh, Hyy/qp->num_cand_nbh), V, D);

		double a = 1.6651*sqrt(D[0][0]);
		double b = 1.6651*sqrt(D[1][1]);

		double rax = V[0][0];
		double ray = V[0][1];

		double phi = atan2(ray,rax);

		unsigned int nbh_iter = 0;
		for (unsigned int j = 0; j < qp->num_cand_nbh; j++) {
			unsigned int jdx = qp->cand_nbh[j];

			double qx = q[jdx]->x;
			double qy = q[jdx]->y;

			double xc = qx - px;
			double yc = qy - py;

			double left  = ((cos(phi)*xc + sin(phi)*yc)/a);
			double right = ((sin(phi)*xc - cos(phi)*yc)/b);

			bool in = left*left + right*right <= 1;

			if (in) {
				qp->nbh[nbh_iter] = jdx;
				nbh_iter++;
			}
		}

		qp->num_nbh = nbh_iter;
		assert(nbh_iter > 0);
	}
}

template<typename T>
void nbh_ellipse(T **p, unsigned int np) {
	nbh_ellipse(p, np, p);
}

#endif /* NEIGHBOR_ELLIPSE_H_ */
