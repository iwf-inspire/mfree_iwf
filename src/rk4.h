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

//runge kutta 4 time integration scheme

#ifndef RK4_H_
#define RK4_H_

#include <stdio.h>
#include <type_traits>

#include "particle.h"
#include "body.h"
#include "simulation_time.h"
#include "correctors.h"
#include "precomp_shape_functions.h"
#include "neighbor_ellipse.h"

//template <typename T> class body;

template<typename T>
class rk4 {
	template<typename U>
	friend class body;
	static_assert(std::is_base_of<particle, T>::value, "T ist not derived from particle");

private:
	T **m_k1 = 0;
	T **m_k2 = 0;
	T **m_k3 = 0;
	T **m_k4 = 0;

	T **m_k1_q = 0;
	T **m_k2_q = 0;
	T **m_k3_q = 0;
	T **m_k4_q = 0;

	template <typename Q=T>
	void rk_add(Q **k_in, Q** p_init, Q** k_out, unsigned int n, double dt) const {
		for (unsigned int i = 0; i < n; i++) {
			p_init[i]->copy_into(k_out[i]);

			k_out[i]->x    += k_in[i]->x_t*dt;
			k_out[i]->y    += k_in[i]->y_t*dt;

			k_out[i]->vx   += k_in[i]->vx_t*dt;
			k_out[i]->vy   += k_in[i]->vy_t*dt;

			k_out[i]->rho  += k_in[i]->rho_t*dt;

			k_out[i]->Sxx  += k_in[i]->Sxx_t*dt;
			k_out[i]->Sxy  += k_in[i]->Sxy_t*dt;
			k_out[i]->Syy  += k_in[i]->Syy_t*dt;
		}
	}

	template <typename Q>
	void step_strong(body<Q> *body) {
		simulation_time *time = &simulation_time::getInstance();
		double m_dt = time->m_dt;

		for (unsigned int i = 0; i < body->m_num_part; i++) {
			body->m_particles[i]->copy_into(m_k1[i]);
		}

		body->compute_time_derivatives(m_k1, body->m_num_part);

		rk_add(m_k1, body->m_particles, m_k2, body->m_num_part, 0.5*m_dt);

		body->compute_time_derivatives(m_k2, body->m_num_part);

		rk_add(m_k2, body->m_particles, m_k3, body->m_num_part, 0.5*m_dt);

		body->compute_time_derivatives(m_k3, body->m_num_part);

		rk_add(m_k3, body->m_particles, m_k4, body->m_num_part, m_dt);

		body->compute_time_derivatives(m_k4, body->m_num_part);

		for (unsigned int i = 0; i < body->m_num_part; i++) {
			body->m_particles[i]->x   += m_dt*1./6.*(m_k1[i]->x_t   + 2*m_k2[i]->x_t   + 2*m_k3[i]->x_t   + m_k4[i]->x_t);
			body->m_particles[i]->y   += m_dt*1./6.*(m_k1[i]->y_t   + 2*m_k2[i]->y_t   + 2*m_k3[i]->y_t   + m_k4[i]->y_t);
			body->m_particles[i]->rho += m_dt*1./6.*(m_k1[i]->rho_t + 2*m_k2[i]->rho_t + 2*m_k3[i]->rho_t + m_k4[i]->rho_t);
			body->m_particles[i]->vx  += m_dt*1./6.*(m_k1[i]->vx_t  + 2*m_k2[i]->vx_t  + 2*m_k3[i]->vx_t  + m_k4[i]->vx_t);
			body->m_particles[i]->vy  += m_dt*1./6.*(m_k1[i]->vy_t  + 2*m_k2[i]->vy_t  + 2*m_k3[i]->vy_t  + m_k4[i]->vy_t);
			body->m_particles[i]->Sxx += m_dt*1./6.*(m_k1[i]->Sxx_t + 2*m_k2[i]->Sxx_t + 2*m_k3[i]->Sxx_t + m_k4[i]->Sxx_t);
			body->m_particles[i]->Sxy += m_dt*1./6.*(m_k1[i]->Sxy_t + 2*m_k2[i]->Sxy_t + 2*m_k3[i]->Sxy_t + m_k4[i]->Sxy_t);
			body->m_particles[i]->Syy += m_dt*1./6.*(m_k1[i]->Syy_t + 2*m_k2[i]->Syy_t + 2*m_k3[i]->Syy_t + m_k4[i]->Syy_t);
		}

		body->m_bc.carry_out_boundary_conditions(body->m_particles, body->m_num_part);

		time->increment_time();
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_ul>::value>::type
	step(body<T> *body) {
		step_strong(body);
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_monaghan>::value>::type
	step(body<T> *body) {
		step_strong(body);
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_tl>::value>::type
	step(body<T> *body) {
		step_strong(body);
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_corotational>::value>::type
	step(body<T> *body) {
		step_strong(body);
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_potential>::value>::type
	step(body<T> *body) {
		step_strong(body);
	}

	template <typename Q = T>
	void do_step_weak(body<T> *body) {
		simulation_time *time = &simulation_time::getInstance();
		double m_dt = time->m_dt;

		for (unsigned int i = 0; i < body->m_num_part; i++) {
			body->m_particles[i]->copy_into(m_k1[i]);
		}

		for (unsigned int i = 0; i < body->m_num_quad; i++) {
			body->m_quad_points[i]->copy_into(m_k1_q[i]);
		}

		body->compute_time_derivatives(m_k1, body->m_num_part, m_k1_q, body->m_num_quad);
		rk_add(m_k1,   body->m_particles,   m_k2,   body->m_num_part, 0.5*m_dt);
		rk_add(m_k1_q, body->m_quad_points, m_k2_q, body->m_num_quad, 0.5*m_dt);

		body->compute_time_derivatives(m_k2, body->m_num_part, m_k2_q, body->m_num_quad);
		rk_add(m_k2,   body->m_particles,   m_k3,   body->m_num_part, 0.5*m_dt);
		rk_add(m_k2_q, body->m_quad_points, m_k3_q, body->m_num_quad, 0.5*m_dt);

		body->compute_time_derivatives(m_k3, body->m_num_part, m_k3_q, body->m_num_quad);
		rk_add(m_k3,   body->m_particles,   m_k4,   body->m_num_part, m_dt);
		rk_add(m_k3_q, body->m_quad_points, m_k4_q, body->m_num_quad, m_dt);

		body->compute_time_derivatives(m_k4, body->m_num_part, m_k4_q, body->m_num_quad);

		for (unsigned int i = 0; i < body->m_num_part; i++) {
			body->m_particles[i]->x   += m_dt*1./6.*(m_k1[i]->x_t   + 2*m_k2[i]->x_t   + 2*m_k3[i]->x_t   + m_k4[i]->x_t);
			body->m_particles[i]->y   += m_dt*1./6.*(m_k1[i]->y_t   + 2*m_k2[i]->y_t   + 2*m_k3[i]->y_t   + m_k4[i]->y_t);
			body->m_particles[i]->rho += m_dt*1./6.*(m_k1[i]->rho_t + 2*m_k2[i]->rho_t + 2*m_k3[i]->rho_t + m_k4[i]->rho_t);

			body->m_particles[i]->vx  += m_dt*1./6.*(m_k1[i]->vx_t  + 2*m_k2[i]->vx_t  + 2*m_k3[i]->vx_t  + m_k4[i]->vx_t);
			body->m_particles[i]->vy  += m_dt*1./6.*(m_k1[i]->vy_t  + 2*m_k2[i]->vy_t  + 2*m_k3[i]->vy_t  + m_k4[i]->vy_t);
		}

		for (unsigned int i = 0; i < body->m_num_quad; i++) {
			body->m_quad_points[i]->x   += m_dt*1./6.*(m_k1_q[i]->x_t   + 2*m_k2_q[i]->x_t   + 2*m_k3_q[i]->x_t   + m_k4_q[i]->x_t);
			body->m_quad_points[i]->y   += m_dt*1./6.*(m_k1_q[i]->y_t   + 2*m_k2_q[i]->y_t   + 2*m_k3_q[i]->y_t   + m_k4_q[i]->y_t);
			body->m_quad_points[i]->rho += m_dt*1./6.*(m_k1_q[i]->rho_t + 2*m_k2_q[i]->rho_t + 2*m_k3_q[i]->rho_t + m_k4_q[i]->rho_t);

			body->m_quad_points[i]->Sxx += m_dt*1./6.*(m_k1_q[i]->Sxx_t + 2*m_k2_q[i]->Sxx_t + 2*m_k3_q[i]->Sxx_t + m_k4_q[i]->Sxx_t);
			body->m_quad_points[i]->Sxy += m_dt*1./6.*(m_k1_q[i]->Sxy_t + 2*m_k2_q[i]->Sxy_t + 2*m_k3_q[i]->Sxy_t + m_k4_q[i]->Sxy_t);
			body->m_quad_points[i]->Syy += m_dt*1./6.*(m_k1_q[i]->Syy_t + 2*m_k2_q[i]->Syy_t + 2*m_k3_q[i]->Syy_t + m_k4_q[i]->Syy_t);
		}
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_ul_weak>::value>::type
	step(body<T> *body) {
		simulation_time *time = &simulation_time::getInstance();

		do_step_weak(body);

		time->increment_time();

		body->m_bc.carry_out_boundary_conditions(body->m_particles, body->m_num_part);

		if (body->get_sim_data().get_correction_constants().do_smooth()) {
			nbh_ellipse(body->m_particles, body->m_num_part);
			body->m_cur_particles = body->m_particles;
			riemann_update_integration_weights(*body);
			precomp_rkpm(body->m_particles, body->m_num_part);
			correctors_smooth_vel<Q>(*body);
		}
	}

	template <typename Q = T>
	typename std::enable_if<std::is_same<Q, particle_tl_weak>::value>::type
	step(body<T> *body) {
		simulation_time *time = &simulation_time::getInstance();

		do_step_weak(body);

		time->increment_time();

		body->m_bc.carry_out_boundary_conditions(body->m_particles, body->m_num_part);
		correctors_smooth_vel<Q>(*body);
	}

public:

	rk4() {}

	rk4(unsigned int num_part, double dt) {
		m_k1 = new T*[num_part];
		m_k2 = new T*[num_part];
		m_k3 = new T*[num_part];
		m_k4 = new T*[num_part];

		for (unsigned int i = 0; i < num_part; i++) {
			m_k1[i] = new T();
			m_k2[i] = new T();
			m_k3[i] = new T();
			m_k4[i] = new T();
		}

		m_k1_q = 0;
		m_k2_q = 0;
		m_k3_q = 0;
		m_k4_q = 0;

		simulation_time *time = &simulation_time::getInstance();
		time->m_dt = dt;
	}

	rk4(unsigned int num_part, unsigned int num_quad, double dt) {
		m_k1 = new T*[num_part];
		m_k2 = new T*[num_part];
		m_k3 = new T*[num_part];
		m_k4 = new T*[num_part];

		for (unsigned int i = 0; i < num_part; i++) {
			m_k1[i] = new T();
			m_k2[i] = new T();
			m_k3[i] = new T();
			m_k4[i] = new T();
		}

		m_k1_q = new T*[num_quad];
		m_k2_q = new T*[num_quad];
		m_k3_q = new T*[num_quad];
		m_k4_q = new T*[num_quad];

		for (unsigned int i = 0; i < num_quad; i++) {
			m_k1_q[i] = new T();
			m_k2_q[i] = new T();
			m_k3_q[i] = new T();
			m_k4_q[i] = new T();
		}

		simulation_time *time = &simulation_time::getInstance();
		time->m_dt = dt;
	}
};

#endif /* RK4_H_ */
