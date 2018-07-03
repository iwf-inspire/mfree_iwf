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

// continuity, momentum and advection equations

#ifndef CONT_MECH_H_
#define CONT_MECH_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
void contmech_momentum_ul(body<T> &b) {
	static_assert(std::is_base_of<particle_ul, T>::value, "T ist not derived from particle_ul");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		double Sxx_x = p->Sxx_x;
		double Sxy_y = p->Sxy_y;
		double Sxy_x = p->Sxy_x;
		double Syy_y = p->Syy_y;

		double rho = p->rho;

		p->vx_t += 1./rho*(Sxx_x + Sxy_y) + p->fcx / p->m;
		p->vy_t += 1./rho*(Sxy_x + Syy_y) + p->fcy / p->m;
	}
}

template<typename T>
void contmech_momentum_ul_weak(body<T> &b) {
	static_assert(std::is_base_of<particle_ul, T>::value, "T ist not derived from particle_ul");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		double Sxx_x = p->Sxx_x;
		double Sxy_y = p->Sxy_y;
		double Sxy_x = p->Sxy_x;
		double Syy_y = p->Syy_y;

		double m = p->m;

		p->vx_t += 1./m*(Sxx_x + Sxy_y) + p->fcx / p->m;
		p->vy_t += 1./m*(Sxy_x + Syy_y) + p->fcy / p->m;
	}
}

template<typename T>
void contmech_momentum_tl(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tt4_ref = load('mfree_iwf/out4.txt');t4_ref = load('mfree_iwf/out4.txt');l");

	double rho0 = b.get_sim_data().get_physical_constants().rho0();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		double Pxx_x = p->Pxx_x;
		double Pxy_x = p->Pxy_x;
		double Pyx_y = p->Pyx_y;
		double Pyy_y = p->Pyy_y;

		p->vx_t += 1./rho0*(Pxx_x + Pyx_y) + p->fcx / p->m;
		p->vy_t += 1./rho0*(Pxy_x + Pyy_y) + p->fcy / p->m;
	}
}

template<typename T>
void contmech_momentum_tl_weak(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tt4_ref = load('mfree_iwf/out4.txt');t4_ref = load('mfree_iwf/out4.txt');l");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		double Pxx_x = p->Pxx_x;
		double Pxy_x = p->Pxy_x;
		double Pyx_y = p->Pyx_y;
		double Pyy_y = p->Pyy_y;

		double m = p->m;

		p->vx_t += 1./m*(Pxx_x + Pyx_y) + p->fcx / p->m;
		p->vy_t += 1./m*(Pxy_x + Pyy_y) + p->fcy / p->m;
	}
}

template<typename T>
void do_contmech_continuity(T **particles, unsigned int n) {

	for (unsigned int i = 0; i < n; i++) {
		T *p = particles[i];

		double rho = p->rho;
		double vx_x = p->vx_x;
		double vy_y = p->vy_y;

		p->rho_t -= rho*(vx_x + vy_y);
	}
}

template<typename T>
void contmech_continuity(body<T> &b) {

	do_contmech_continuity(b.get_cur_particles(), b.get_num_part());
}

template<typename T>
void contmech_continuity_gp(body<T> &b) {

	do_contmech_continuity(b.get_cur_quad_points(), b.get_num_quad_points());
}

template<typename T>
void contmech_advection(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];

		p->x_t += p->vx;
		p->y_t += p->vy;
	}
}

template<typename T>
void contmech_advection_gp(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *p = b.get_cur_quad_points()[i];

		p->x_t += p->vx;
		p->y_t += p->vy;
	}
}

template<typename T>
void do_contmech_cauchy_to_nominal(T **particles, unsigned int n) {
	for (unsigned int i = 0; i < n; i++) {
		T *p = particles[i];

		glm::dmat2x2 F(p->Hxx+1, p->Hyx, p->Hxy, p->Hyy+1);
		double J = glm::determinant(F);
		glm::dmat2x2 sigma = glm::dmat2x2(p->Sxx - p-> p, p->Sxy, p->Sxy, p->Syy - p->p);
		glm::dmat2x2 Finv = glm::inverse(F);
		glm::dmat2x2 P = J*Finv*sigma;
		p->Pxx = P[0][0];
		p->Pxy = P[1][0];
		p->Pyx = P[0][1];
		p->Pyy = P[1][1];
	}
}

template<typename T>
void contmech_velocity_gradient_tl_gp(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");
	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *p = b.get_cur_quad_points()[i];
		glm::dmat2x2 F(p->Hxx+1., p->Hyx, p->Hxy, p->Hyy+1.);
		glm::dmat2x2 Fdot(p->Fdotxx, p->Fdotyx, p->Fdotxy, p->Fdotyy);
		glm::dmat2x2 Finv = glm::inverse(F);
		glm::dmat2x2 L = Fdot*Finv;

		p->vx_x = L[0][0];
		p->vx_y = L[1][0];
		p->vy_x = L[0][1];
		p->vy_y = L[1][1];
	}
}

template<typename T>
void contmech_velocity_gradient_tl(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *p = b.get_cur_particles()[i];
		glm::dmat2x2 F(p->Hxx+1., p->Hyx, p->Hxy, p->Hyy+1.);
		glm::dmat2x2 Fdot(p->Fdotxx, p->Fdotyx, p->Fdotxy, p->Fdotyy);
		glm::dmat2x2 Finv = glm::inverse(F);
		glm::dmat2x2 L = Fdot*Finv;

		p->vx_x = L[0][0];
		p->vx_y = L[1][0];
		p->vy_x = L[0][1];
		p->vy_y = L[1][1];
	}
}


template<typename T>
void contmech_cauchy_to_nominal_gp(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	do_contmech_cauchy_to_nominal(b.get_cur_quad_points(), b.get_num_quad_points());
}

template<typename T>
void contmech_cauchy_to_nominal(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	do_contmech_cauchy_to_nominal(b.get_cur_particles(), b.get_num_part());
}

#endif /* CONT_MECH_H_ */
