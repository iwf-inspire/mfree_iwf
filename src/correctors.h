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

//compute artificial stress, artificial viscosity, XSPH
// all for example given in Gray, J. P., J. J. Monaghan, and R. P. Swift. "SPH elastic dynamics." Computer methods in applied mechanics and engineering 190.49-50 (2001): 6641-6662.
//
// velocity field smoothing according to
//Ostad-Hossein, Hassan, and Soheil Mohammadi. "A field smoothing stabilization of particle methods in elastodynamics." Finite Elements in Analysis and Design 44.9-10 (2008): 564-579.

#ifndef CORRECTORS_H_
#define CORRECTORS_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"
#include "body.h"

template <typename T> class body;

static double stress_angle(double sxx, double sxy, double syy, double eps) {
	double numer = 2.*sxy;
	double denom = sxx - syy + eps;
	return 0.5*atan2(numer,denom);
}

template<typename T>
void correctors_mghn_artificial_stress(body<T> &b) {
	static_assert(std::is_base_of<particle_monaghan, T>::value, "T ist not derived from particle_monaghan");

	double eps = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_eps();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		double rhoi = p->rho;

		// diag
		double sxx = p->Sxx - p->p;
		double syy = p->Syy - p->p;

		// off diag
		double sxy = p->Sxy;

		double theta = stress_angle(sxx,sxy,syy,0.);

		double cos_theta = cos(theta);
		double sin_theta = sin(theta);

		double cos_theta2 = cos_theta*cos_theta;
		double sin_theta2 = sin_theta*sin_theta;

		double rot_sxx = cos_theta2*sxx + 2.0*cos_theta*sin_theta*sxy + sin_theta2*syy;
		double rot_syy = sin_theta2*sxx - 2.0*cos_theta*sin_theta*sxy + cos_theta2*syy;

		double rot_rxx = 0.;
		double rot_ryy = 0.;
		if (rot_sxx > 0) rot_rxx = -eps*rot_sxx/(rhoi*rhoi);
		if (rot_syy > 0) rot_ryy = -eps*rot_syy/(rhoi*rhoi);

		p->Rxx = cos_theta2*rot_rxx + sin_theta2*rot_ryy;
		p->Rxy = cos_theta*sin_theta*(rot_rxx - rot_ryy);
		p->Ryy = sin_theta2*rot_rxx + cos_theta2*rot_ryy;
	}
}

template<typename T>
void correctors_mghn_artificial_viscosity(body<T> &b) {
	double alpha = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_alpha();
	double beta  = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_beta();
	double eta   = b.get_sim_data().get_correction_constants().get_art_visc_const().artvisc_eta();

	double K = b.get_sim_data().get_physical_constants().K();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double xi   = pi->x;
		double yi   = pi->y;
		double vxi  = pi->vx;
		double vyi  = pi->vy;
		double hi   = pi->h;
		double rhoi = pi->rho;
		double ci   = sqrt(K/rhoi);

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			particle *pj = b.get_cur_particles()[jdx];

			double xj   = pj->x;
			double yj   = pj->y;
			double vxj  = pj->vx;
			double vyj  = pj->vy;
			double hj   = pj->h;
			double mj   = pj->m;
			double rhoj = pj->rho;
			double cj   = sqrt(K/rhoj);

			double vijx = vxi - vxj;
			double vijy = vyi - vyj;

			double xij = xi - xj;
			double yij = yi - yj;

			double vijposij = vijx*xij + vijy*yij;

			if (vijposij >= 0.) continue;

			double cij = 0.5*(ci+cj);
			double hij = 0.5*(hi+hj);
			double rhoij = 0.5*(rhoi+rhoj);
			double r2ij = xij*xij + yij*yij;
			double muij = (hij*vijposij)/(r2ij + eta*eta*hij*hij);
			double piij = (-alpha*cij*muij + beta*muij*muij)/rhoij;

			pi->vx_t += -mj*piij*w.w_x;
			pi->vy_t += -mj*piij*w.w_y;
		}
	}
}

template<typename T>
void correctors_xsph(body<T> &b) {
	double eps = b.get_sim_data().get_correction_constants().xsph_eps();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double vxi = pi->vx;
		double vyi = pi->vy;
		double rhoi = pi->rho;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj  = pj->vx;
			double vyj  = pj->vy;
			double rhoj = pj->rho;

			double vijx = vxi - vxj;
			double vijy = vyi - vyj;

			double rhoij = 0.5*(rhoi + rhoj);
			double fac = -eps*w.w*pi->m/rhoij;

			pi->x_t += fac*vijx;
			pi->y_t += fac*vijy;
		}
	}
}

template<typename T>
void correctors_smooth_vel(body<T> &b) {
	if (!b.get_sim_data().get_correction_constants().do_smooth()) return;

	static double *vsx = 0;
	static double *vsy = 0;

	if (vsx == 0) {
		vsx = (double *) calloc(b.get_num_part(), sizeof(double));
		vsy = (double *) calloc(b.get_num_part(), sizeof(double));
	}

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_particles()[i];

		double vsxi = 0.;
		double vsyi = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_particles()[jdx];

			double vxj  = pj->vx;
			double vyj  = pj->vy;
			double quad_weight = pj->quad_weight;

			vsxi += vxj*w.w*quad_weight;
			vsyi += vyj*w.w*quad_weight;
		}

		vsx[i] = vsxi;
		vsy[i] = vsyi;
	}

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		b.get_particles()[i]->vx = vsx[i];
		b.get_particles()[i]->vy = vsy[i];
	}
}

#endif /* CORRECTORS_H_ */
