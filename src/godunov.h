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

//modified momentum and continuity equation for godunov type SPH
// acoustic approximations according to Parshikov, Anatoly N., and Stanislav A. Medin. "Smoothed particle hydrodynamics using interparticle contact algorithms." Journal of computational physics 180.1 (2002): 358-382.

#ifndef GODUNOV_H_
#define GODUNOV_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

struct trafo {
	glm::dmat2x2 l;
	glm::dvec2 sigmaistarR;
	glm::dvec2 sigmajstarR;
	glm::dvec2 vistar;
	glm::dvec2 vjstar;
};

static void riemann_solver_acoustic(double rhol, double rhor,
		double pl, double pr,
		double ul, double ur,
		double *pstar, double *ustar, simulation_data data) {

	double K = data.get_physical_constants().K();

	// Lagrangian sound speeds
	double csl = sqrt(K/rhol);
	double csr = sqrt(K/rhor);

	double frac = 1./(rhol*csl + rhor*csr);

	*ustar = frac*(ul*rhol*csl + ur*rhor*csr - pl + pr);
	*pstar = frac*(pr*rhol*csl + pl*rhor*csr - rhol*csl*rhor*csr*(ul-ur));
}

static trafo godunov_get_trafo(const particle_ul *pi, const particle_ul *pj) {
	//using notation from appendix B, parshikov paper
	double Rx = pj->x - pi->x;
	double Ry = pj->y - pi->y;

	double Sx =  Ry;
	double Sy = -Rx;

	glm::dvec2 R = glm::normalize(glm::dvec2(Rx, Ry));
	glm::dvec2 S = glm::normalize(glm::dvec2(Sx, Sy));

	glm::dmat2x2 sigmai(pi->Sxx - pi->p, pi->Sxy, pi->Sxy, pi->Syy);
	glm::dmat2x2 sigmaj(pj->Sxx - pj->p, pj->Sxy, pj->Sxy, pj->Syy);

	trafo t;
	t.l = glm::mat2x2(R,S);

	glm::dvec2 sigmaiR = sigmai*R;	// this is nothing but carrying out part of the derivative (multiplying by rij/|rij|
	glm::dvec2 sigmajR = sigmaj*R;  // which is part of Wij'rij/(h*|rji|), which can also be seen as projection sigma onto R

	t.sigmaistarR = t.l*sigmaiR;	// this is the actual trafo into system RS
	t.sigmajstarR = t.l*sigmajR;

	glm::dvec2 vi(pi->vx, pi->vy);
	glm::dvec2 vj(pj->vx, pj->vy);

	t.vistar = t.l*vi;
	t.vjstar = t.l*vj;

	return t;
}

template<typename T>
void godunov_momentum(body<T> &b) {
	static_assert(std::is_base_of<particle_ul, T>::value, "T ist not derived from particle_ul");

	// variables needed for riemann solver
	double K = b.get_sim_data().get_physical_constants().K();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double ci = sqrt(K/pi->rho);
		double rhoi = pi->rho;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			if (pi->idx == pj->idx) continue;

			trafo t = godunov_get_trafo(pi, pj);

			double cj = sqrt(K/pj->rho);
			double rhoj = pj->rho;
			double quad_weight = pj->quad_weight;

			// Riemann solver (Parshikov Paper, (23)-(28))
			glm::dvec2 SRijstar;
			double Riemann_frac = 1./(rhoi*ci + rhoj*cj);
			SRijstar[0] = Riemann_frac*(t.sigmajstarR[0]*rhoi*ci + t.sigmaistarR[0]*rhoj*cj + rhoi*ci*rhoj*cj*(t.vjstar[0] - t.vistar[0]));
			SRijstar[1] = Riemann_frac*(t.sigmajstarR[1]*rhoi*ci + t.sigmajstarR[1]*rhoj*cj + rhoi*ci*rhoj*cj*(t.vjstar[1] - t.vistar[1]));

			//parshikov paper, appendix b, l should be transposed imho
			glm::dvec2 SRij = glm::transpose(t.l)*SRijstar;

			double rx = pj->x - pi->x;
			double ry = pj->y - pi->y;
			double mag_r = sqrt(rx*rx + ry*ry);

			if (fabs(rx) > 0.) {
				pi->vx_t += 2./rhoi*w.w_x*mag_r/rx*SRij[0]*quad_weight;
			}

			if (fabs(ry) > 0.) {
				pi->vy_t += 2./rhoi*w.w_y*mag_r/ry*SRij[1]*quad_weight;
			}
		}
	}

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];
		pi->vx_t += pi->fcx/pi->m;
		pi->vy_t += pi->fcy/pi->m;
	}
}

template<typename T>
void godunov_continuity(body<T> &b) {
	static_assert(std::is_base_of<particle_ul, T>::value, "T ist not derived from particle_ul");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double xi   = pi->x;
		double yi   = pi->y;
		double rhoi = pi->rho;

		double rho_t = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			if (pi->idx == pj->idx) continue;

			double xj = pj->x;
			double yj = pj->y;

			double quad_weight = pj->quad_weight;

			double dx = xi - xj;
			double dy = yi - yj;
			double mag_r = sqrt(dx*dx + dy*dy);

			double nx = dx/mag_r;
			double ny = dy/mag_r;

			double vir = pi->vx*nx + pi->vy*ny;
			double vjr = pj->vx*nx + pj->vy*ny;

			double pstar, vstar;
			riemann_solver_acoustic(pi->rho, pj->rho, pi->p, pj->p, vir, vjr, &pstar, &vstar, b.get_sim_data());

			rho_t += 2*rhoi*(vir - vstar)*(w.w_x*nx + w.w_y*ny)*quad_weight;
		}

		pi->rho_t = rho_t;
	}
}

#endif /* GODUNOV_H_ */
