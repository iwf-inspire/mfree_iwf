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

#ifndef BODY_H_
#define BODY_H_

#include <type_traits>
#include <stdio.h>
#include <vector>
#include <algorithm>

#include "grid.h"
#include "particle.h"
#include "simulation_data.h"
#include "precomp_shape_functions.h"
#include "rk4.h"
#include "boundary_conditions.h"
#include "boundary.h"


template <typename T>
class body {
	template<typename U>
	friend class rk4;
	static_assert(std::is_base_of<particle, T>::value, "T ist not derived from particle");

private:
	grid<T> m_grid;
	simulation_data m_data;
	std::vector< void (*)(body<T> &b) > m_actions;
	boundary_conditions<T> m_bc;
	boundary<T> *m_boundary = 0;

	unsigned int m_num_part;
	T** m_particles;
	T** m_cur_particles;

	unsigned int m_num_quad;
	T** m_quad_points;
	T** m_cur_quad_points;

	rk4<T> m_time_stepper;

	gauss_int *m_gauss_integration = 0;

	template <typename Q = T>
	typename std::enable_if<std::is_base_of<particle_tl_weak, Q>::value>::type
	compute_time_derivatives(T** particles, unsigned int num_part, T** quad_points, unsigned int num_quad) {
		for (unsigned int i = 0; i < num_part; i++) {
			particles[i]->reset();
		}

		for (unsigned int i = 0; i < num_quad; i++) {
			quad_points[i]->reset();
		}

		m_bc.carry_out_boundary_conditions(particles, num_part);

		m_cur_particles = particles;
		m_cur_quad_points = quad_points;
		for (auto it = m_actions.begin(); it != m_actions.end(); ++it) {
			(*it)(*this);
		}
	}

	template <typename Q = T>
	typename std::enable_if<std::is_base_of<particle_ul_weak, Q>::value>::type
	compute_time_derivatives(T** particles, unsigned int num_part, T** quad_points, unsigned int num_quad) {
		for (unsigned int i = 0; i < num_part; i++) {
			particles[i]->reset();
		}

		for (unsigned int i = 0; i < num_quad; i++) {
			quad_points[i]->reset();
		}

		if (!m_gauss_integration) {
			nbh_ellipse(particles, m_num_part);
			nbh_ellipse(quad_points, m_num_quad, particles);

			precomp_rkpm<Q>(particles, quad_points, num_quad);
			precomp_rkpm<Q>(particles, num_part);
		} else {
			m_gauss_integration->update_gauss_points(particles, quad_points);
		}

		m_cur_particles = particles;
		m_cur_quad_points = quad_points;
		for (auto it = m_actions.begin(); it != m_actions.end(); ++it) {
			(*it)(*this);
		}

		m_bc.carry_out_boundary_conditions(particles, num_part);
	}

	template <typename Q = T>
	typename std::enable_if<std::is_base_of<particle_tl, Q>::value>::type
	compute_time_derivatives(T** particles, unsigned int num_part) {
		for (unsigned int i = 0; i < num_part; i++) {
			particles[i]->reset();
		}

		m_bc.carry_out_boundary_conditions(particles, num_part);

		m_cur_particles = particles;
		for (auto it = m_actions.begin(); it != m_actions.end(); ++it) {
			(*it)(*this);
		}
	}

	template <typename Q = T>
	typename std::enable_if<std::is_base_of<particle_ul, Q>::value>::type
	compute_time_derivatives(T** particles, unsigned int num_part) {
		for (unsigned int i = 0; i < num_part; i++) {
			particles[i]->reset();
		}

//		m_bc.carry_out_boundary_conditions(particles, num_part);		//BAD IDEA!

		m_grid.update_geometry(particles, num_part, 2.);
		m_grid.assign_hashes(particles, num_part);

		std::sort(particles, particles+num_part, [](const T *a, const T *b) {return a->hash < b->hash;});
		m_grid.construct_verlet_lists(particles, num_part, 2.);

		precomp_sph(particles, num_part);

		m_cur_particles = particles;
		for (auto it = m_actions.begin(); it != m_actions.end(); ++it) {
			(*it)(*this);
		}

		std::sort(particles, particles+num_part, [](const T *a, const T *b) {return a->idx < b->idx;});
	}


public:

	void update_boundary() {
		if (!m_boundary) return;

		kd_point *boundary = m_boundary->points();
		const unsigned int *indexes = m_boundary->indexes();

		for (unsigned int i = 0; i < m_boundary->num_bound(); i++) {
			unsigned int idx = indexes[i];
			boundary[i].x = m_particles[idx]->x;
			boundary[i].y = m_particles[idx]->y;
		}

		m_boundary->update_normals();
	}

	void debug_dump_boundary(unsigned int identifier, unsigned int step) {
		char buf[256];
		sprintf(buf, "boundary_%04d_%04d.txt", step, identifier);
		FILE *fp = fopen(buf, "w+");
		for (unsigned int i = 0; i < m_boundary->num_bound(); i++) {
			fprintf(fp, "%f %f %f %f\n", m_boundary->points()[i].x, m_boundary->points()[i].y, m_boundary->points()[i].nx, m_boundary->points()[i].ny);
		}
		fclose(fp);
	}

	void set_boundary(boundary<T> *boundary) {
		m_boundary = boundary;
	}

	const boundary<T> *get_boundary() const {
		return (const boundary<T>*) m_boundary;
	}

	bool has_boundary() const {
		return m_boundary != NULL;
	}

	simulation_data get_sim_data() const {
		return m_data;
	}

	T** get_cur_particles() const {
		return m_cur_particles;
	}

	T** get_particles() const {
		return m_particles;
	}

	unsigned int get_num_part () const {
		return m_num_part;
	}

	T** get_cur_quad_points() const {
		return m_cur_quad_points;
	}

	unsigned int get_num_quad_points () const {
		return m_num_quad;
	}

	body() {};

	body(T** particles, unsigned int n, simulation_data data, double dt, boundary_conditions<T> bc = boundary_conditions<T>(),
			T** quad_points = 0, unsigned int nq = 0, gauss_int *gauss_integartion = 0);

	void add_action(void (*action) (body<T> &b)) {
		m_actions.push_back(action);
	}

	void step() {
		m_time_stepper.step(this);
	}
};

template<typename T>
body<T>::body(T** particles, unsigned int n, simulation_data data, double dt, boundary_conditions<T> bc, T** quad_points, unsigned int nq, gauss_int *gauss_integartion) :
	m_data(data), m_bc(bc),
	m_num_part(n), m_particles(particles), m_cur_particles(0),
	m_num_quad(nq), m_quad_points(quad_points), m_cur_quad_points(0),
	m_time_stepper(n, nq, dt), m_gauss_integration(gauss_integartion) {}

#endif /* BODY_H_ */
