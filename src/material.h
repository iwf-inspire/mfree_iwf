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

// material modelling
//	standard eos + jaumann rate
//	small strain "state" form (cauchy green)
//	large strain "rate" form (green st venant)
//	extremely simple perfect plasticty as given in Gray, J. P., J. J. Monaghan, and R. P. Swift. "SPH elastic dynamics." Computer methods in applied mechanics and engineering 190.49-50 (2001): 6641-6662.


#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
static void do_material_eos(T** particles, unsigned int n, double rho0, double K) {
	for (unsigned int i = 0; i < n; i++) {
		T *p = particles[i];
		double c0 = sqrt(K/rho0);
		p->p = c0*c0*(p->rho - rho0);
	}
}

template<typename T>
void material_eos_gp(body<T> &b) {
//	static_assert(std::is_base_of<particle_ul_weak, T>::value, "T ist not derived from particle_ul");

	double rho0 = b.get_sim_data().get_physical_constants().rho0();
	double K    = b.get_sim_data().get_physical_constants().K();

	do_material_eos(b.get_cur_quad_points(), b.get_num_quad_points(), rho0, K);
}

template<typename T>
void material_eos(body<T> &b) {
//	static_assert(std::is_base_of<particle_ul, T>::value, "T ist not derived from particle_ul");

	double rho0 = b.get_sim_data().get_physical_constants().rho0();
	double K    = b.get_sim_data().get_physical_constants().K();

	do_material_eos(b.get_cur_particles(), b.get_num_part(), rho0, K);
}

template<typename T>
static void do_material_stress_rate_jaumann(T** particles, unsigned int n, double G) {
	for (unsigned int i = 0; i < n; i++) {
		T *p = particles[i];

		glm::dmat2x2 epsdot = glm::dmat2x2(p->vx_x, 0.5*(p->vx_y + p->vy_x), 0.5*(p->vx_y + p->vy_x), p->vy_y);
		glm::dmat2x2 omega  = glm::dmat2x2(0.     , 0.5*(p->vy_x - p->vx_y), 0.5*(p->vx_y - p->vy_x), 0.     );
		glm::dmat2x2 S      = glm::dmat2x2(p->Sxx, p->Sxy, p->Sxy, p->Syy);
		glm::dmat2x2 I      = glm::dmat2x2(1.);

		double trace_epsdot = epsdot[0][0] + epsdot[1][1];

		glm::dmat2x2 S_t = 2*G*(epsdot - 1./3.*trace_epsdot*I) - S*omega + omega*S;

		p->Sxx_t += S_t[0][0];
		p->Sxy_t += S_t[0][1];
		p->Syy_t += S_t[1][1];
	}
}

template<typename T>
void material_stress_rate_jaumann_gp(body<T> &b) {
//	static_assert(std::is_base_of<particle_ul_weak, T>::value, "T ist not derived from particle_ul");

	double G = b.get_sim_data().get_physical_constants().G();

	do_material_stress_rate_jaumann(b.get_cur_quad_points(), b.get_num_quad_points(), G);
}

template<typename T>
void material_stress_rate_jaumann(body<T> &b) {
	double G = b.get_sim_data().get_physical_constants().G();

	do_material_stress_rate_jaumann(b.get_cur_particles(), b.get_num_part(), G);
}

template<typename T>
void do_material_compute_stress_cauchy_green(T** particles, unsigned int  n, double E, double nu) {
	double D  = E/(1-nu*nu);

	for (unsigned int i = 0; i < n; i++) {
		glm::dmat2x2 F = glm::dmat2x2(particles[i]->Hxx+1., particles[i]->Hyx, particles[i]->Hxy, particles[i]->Hyy+1.);
		glm::dmat2x2 Ecg = 0.5*(F + glm::transpose(F)) - glm::dmat2x2(1.); //cauchy green

		particles[i]->Sxx = D*(Ecg[0][0]    + nu*Ecg[1][1]);
		particles[i]->Syy = D*(Ecg[0][0]*nu +    Ecg[1][1]);
		particles[i]->Sxy = D*(1-nu)*Ecg[0][1];
	}
}

template<typename T>
void material_compute_stress_cauchy_green_gp(body<T> &b) {
	double E  = b.get_sim_data().get_physical_constants().E();
	double nu = b.get_sim_data().get_physical_constants().nu();

	do_material_compute_stress_cauchy_green(b.get_cur_quad_points(), b.get_num_quad_points(), E, nu);
}

template<typename T>
void material_compute_stress_cauchy_green(body<T> &b) {
	double E  = b.get_sim_data().get_physical_constants().E();
	double nu = b.get_sim_data().get_physical_constants().nu();

	do_material_compute_stress_cauchy_green(b.get_cur_particles(), b.get_num_part(), E, nu);
}

template<typename T>
static void do_material_compute_stress_green_stvenant(T** particles, unsigned int  n, double E, double nu) {
	double D  = E/(1-nu*nu);

	for (unsigned int i = 0; i < n; i++) {
		glm::dmat2x2 F = glm::dmat2x2(particles[i]->Hxx+1., particles[i]->Hyx, particles[i]->Hxy, particles[i]->Hyy+1.);
		glm::dmat2x2 C = glm::transpose(F)*F;
		glm::dmat2x2 Esv = 0.5*(C-glm::dmat2x2(1.)); //green st venant

		particles[i]->Sxx = D*(Esv[0][0]    + nu*Esv[1][1]);
		particles[i]->Syy = D*(Esv[0][0]*nu +    Esv[1][1]);
		particles[i]->Sxy = D*(1-nu)*Esv[1][0];
	}
}

template<typename T>
void material_compute_stress_green_stvenant_gp(body<T> &b) {
	double E  = b.get_sim_data().get_physical_constants().E();
	double nu = b.get_sim_data().get_physical_constants().nu();

	do_material_compute_stress_green_stvenant(b.get_cur_quad_points(), b.get_num_quad_points(), E, nu);
}

template<typename T>
void material_compute_stress_green_stvenant(body<T> &b) {
	double E  = b.get_sim_data().get_physical_constants().E();
	double nu = b.get_sim_data().get_physical_constants().nu();

	do_material_compute_stress_green_stvenant(b.get_cur_particles(), b.get_num_part(), E, nu);
}

template<typename T>
void material_plasticity_gray2000(body<T> &b) {
	T **particles = b.get_cur_particles();

	double yvm = 880e6;

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		double sxx = particles[i]->Sxx;
		double sxy = particles[i]->Sxy;
		double syx = particles[i]->Sxy;
		double syy = particles[i]->Syy;

		double J2 = sxx*sxx + syy*syy + sxy*syx + sxx*syy;

		if (J2  < 1e-6) continue;

		double f = fmin(sqrt(((yvm*yvm)/3)/J2),1.);

		particles[i]->Sxx = f*sxx;
		particles[i]->Sxy = f*sxy;
		particles[i]->Syy = f*syy;
	}
}

#endif /* MATERIAL_H_ */
