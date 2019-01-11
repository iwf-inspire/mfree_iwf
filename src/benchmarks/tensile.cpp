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

#include "tensile.h"

template<typename T>
static void left_bndry_fun_tl(T *p) {
	p->x = p->X - 10*simulation_time::getInstance().get_time();
	p->y = p->Y;

	p->vx = -10.;
	p->vy =   0.;

	p->vx_t = 0.;
	p->vy_t = 0.;
}

template<typename T>
static void right_bndry_fun_tl(T *p) {
	p->x = p->X + 10*simulation_time::getInstance().get_time();
	p->y = p->Y;

	p->vx =  10.;
	p->vy =   0.;

	p->vx_t = 0.;
	p->vy_t = 0.;
}

template<typename T>
static void left_bndry_fun_ul(T *p) {
	p->x_t =  -10.;
	p->y_t =   0.;

	p->vx  =  -10.;
	p->vy  =   0.;

	p->vx_t = 0.;
	p->vy_t = 0.;
}

template<typename T>
static void right_bndry_fun_ul(T *p) {
	p->x_t =  +10.;
	p->y_t =   0.;

	p->vx  =  +10.;
	p->vy  =   0.;

	p->vx_t = 0.;
	p->vy_t = 0.;
}

std::vector<body<particle_ul>> tensile_godunov(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1;
	double dx = L/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-6;

	particle_ul **particles  = new particle_ul*[nbox*nbox];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;

			particles[part_iter] = new particle_ul(part_iter);
			particles[part_iter]->x = px;
			particles[part_iter]->y = py;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}

			part_iter++;
		}
	}

	unsigned int n = nbox*nbox;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	boundary_conditions<particle_ul> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_ul);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_ul);

	simulation_data sim_data(physical_constants, correction_constants());

	body<particle_ul> b = body<particle_ul>(particles,  n, sim_data, dt, bc);

	b.add_action(&riemann_update_integration_weights);
	b.add_action(&material_eos);
	b.add_action(&derive_velocity);
	b.add_action(&material_stress_rate_jaumann);
	b.add_action(&godunov_continuity);
	b.add_action(&godunov_momentum);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_ul>>({b});
}

std::vector<body<particle_corotational>> tensile_corot(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);
	double hdx = 1.;

	double dt = 1e-6;

	particle_corotational **particles = new particle_corotational*[nbox*nbox];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;

			particles[part_iter] = new particle_corotational(part_iter);
			particles[part_iter]->x = px;
			particles[part_iter]->y = py;

			particles[part_iter]->X = px;
			particles[part_iter]->Y = py;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}

			part_iter++;

		}
	}

	unsigned int n = nbox*nbox;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	boundary_conditions<particle_corotational> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_tl);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_tl);

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(1., 1., 0.1), 0., false));

	find_neighbors(particles, n, hdx*dx, distance_euclidian);
	precomp_sph<particle_corotational>(particles,  n);

	body<particle_corotational> b = body<particle_corotational>(particles, n, sim_data, dt, bc);

	b.add_action(&corotational_anisotropy_tensor);
	b.add_action(&derive_displacement_corotational);
	b.add_action(&material_compute_stress_cauchy_green);
	b.add_action(&contmech_cauchy_to_nominal);
	b.add_action(&derive_stress_corotational);
	b.add_action(&correctors_mghn_artificial_viscosity);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_corotational>>({b});
}

std::vector<body<particle_tl>> tensile_tl(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);
	double hdx = 1.;

	double dt = 1e-6;

	particle_tl **particles  = new particle_tl*[nbox*nbox];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particles[part_iter] = new particle_tl(part_iter);
			particles[part_iter]->x = px;
			particles[part_iter]->y = py;

			particles[part_iter]->X = px;
			particles[part_iter]->Y = py;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}

			part_iter++;
		}
	}

	unsigned int n = nbox*nbox;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	boundary_conditions<particle_tl> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_tl);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_tl);

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(1., 1., 0.1), 0., false));

	find_neighbors(particles, n, hdx*dx, distance_euclidian);
	precomp_randles_libersky<particle_tl>(particles,  n);

	body<particle_tl> b = body<particle_tl>(particles,  n, sim_data, dt, bc);

	b.add_action(&derive_displacement);
	b.add_action(&derive_velocity_tl);
	b.add_action(&contmech_velocity_gradient_tl);
	b.add_action(&contmech_continuity);
	b.add_action(&material_eos);
	b.add_action(&material_stress_rate_jaumann);
	b.add_action(&contmech_cauchy_to_nominal);
	b.add_action(&derive_stress_tl);
	b.add_action(&correctors_mghn_artificial_viscosity);
	b.add_action(&contmech_momentum_tl);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_tl>>({b});
}

std::vector<body<particle_potential>>    tensile_potential(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);
	double hdx = 1.;

	double dt = 1e-6;

	particle_potential **particles  = new particle_potential*[nbox*nbox];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particles[part_iter] = new particle_potential(part_iter);
			particles[part_iter]->x = px;
			particles[part_iter]->y = py;

			particles[part_iter]->X = px;
			particles[part_iter]->Y = py;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}

			part_iter++;
		}
	}

	unsigned int n = nbox*nbox;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(1., 1., 0.1), 0., false));

	find_neighbors(particles, n, hdx*dx, distance_euclidian);
	precomp_randles_libersky_symmetric<particle_potential>(particles,  n);

	boundary_conditions<particle_potential> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_tl);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_tl);

	body<particle_potential> b = body<particle_potential>(particles,  n, sim_data, dt, bc);

	b.add_action(&derive_displacement);
	b.add_action(&material_compute_stress_green_stvenant);
	b.add_action(&potential_comp_elastic_potential);
	b.add_action(&correctors_mghn_artificial_viscosity);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_potential>>({b});
}

std::vector<body<particle_monaghan>> tensile_monaghan(unsigned int nbox) {

	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);
	double hdx = 1.3;

	double dt = 1e-6;

	particle_monaghan **particles  = new particle_monaghan*[nbox*nbox];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particles[part_iter] = new particle_monaghan(part_iter);
			particles[part_iter]->x = px;
			particles[part_iter]->y = py;

			if (i < 3) {
				left_bndry.push_back(part_iter);
			}

			if (i >= nbox-3) {
				right_bndry.push_back(part_iter);
			}

			part_iter++;
		}
	}

	unsigned int n = nbox*nbox;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	//correction constants
	double alpha = 1.;
	double beta  = 1.;
	double eta   = 0.1;

	double art_stress_eps = 0.3;
	kernel_result w = cubic_spline(0, 0, dx, 0, hdx*dx);
	double wdeltap = w.w;
	double stress_exponent = 4.;

	double xsph_eps = 0.5;

	correction_constants correction_constants(constants_monaghan(wdeltap, stress_exponent, art_stress_eps),
			constants_artificial_viscosity(alpha, beta, eta), xsph_eps);

	simulation_data sim_data(physical_constants, correction_constants);

	boundary_conditions<particle_monaghan> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_ul);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_ul);


	body<particle_monaghan> b = body<particle_monaghan>(particles,  n, sim_data, dt, bc);

	b.add_action(&riemann_update_integration_weights);
	b.add_action(&material_eos);
	b.add_action(&correctors_mghn_artificial_stress);
	b.add_action(&derive_stress_monaghan);
	b.add_action(&derive_velocity);
	b.add_action(&correctors_mghn_artificial_viscosity);
	b.add_action(&correctors_xsph);
	b.add_action(&material_stress_rate_jaumann);
	b.add_action(&contmech_continuity);
	b.add_action(&contmech_momentum_ul);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_monaghan>>({b});
}

std::vector<body<particle_ul_weak>> tensile_ul_fem(unsigned int nbox) {

	// material constants
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);

	double dt = 1e-6;	//stability limit for rk4

	particle_ul_weak **particles = new particle_ul_weak*[nbox*nbox];
	unsigned n = nbox*nbox;

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particle_ul_weak *p = new particle_ul_weak(part_iter);

			p->x = px;
			p->y = py;

			p->X = px;
			p->Y = py;

			particles[part_iter] = p;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}
			part_iter++;
		}
	}

	std::vector<element> elements = element_map_rectangle(nbox-1, nbox-1);
	unsigned num_gp = 2*2*elements.size();
	particle_ul_weak **gauss_points = new particle_ul_weak*[num_gp];
	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i] = new particle_ul_weak(i);
	}
	gauss_int *gi = new gauss_int(elements, 2);
	gi->update_gauss_points(particles, gauss_points);

	boundary_conditions<particle_ul_weak> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_ul);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_ul);

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->quad_weight = 1.;
		particles[i]->m = dx*dx*rho0;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i]->rho = rho0;
	}

	std::vector<unsigned int> boundary_idx = convex_hull(particles, n);
	boundary<particle_ul_weak> *bnd = new boundary<particle_ul_weak>(boundary_idx, particles);

	simulation_data sim_data(physical_constants, correction_constants());
	body<particle_ul_weak> b(particles, n, sim_data, dt, bc, gauss_points, num_gp, gi);

	b.add_action(&material_eos_gp);
	b.add_action(&derive_stress_weak_ul);
	b.add_action(&derive_velocity_weak);
	b.add_action(&material_stress_rate_jaumann_gp);
	b.add_action(&contmech_momentum_ul_weak);
	b.add_action(&contmech_continuity_gp);
	b.add_action(&contmech_advection);

	b.set_boundary(bnd);

	return std::vector<body<particle_ul_weak>>({b});
}

std::vector<body<particle_tl_weak>> tensile_tl_fem(unsigned int nbox) {

	// material constants
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1.;

	//problem dimensions
	double L = 1.;
	double dx = L/(nbox-1);

	double dt = 1e-6;

	particle_tl_weak **particles = new particle_tl_weak*[nbox*nbox];
	unsigned n = nbox*nbox;

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particle_tl_weak *p = new particle_tl_weak(part_iter);

			p->x = px;
			p->y = py;

			p->X = px;
			p->Y = py;

			particles[part_iter] = p;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}
			part_iter++;
		}
	}

	std::vector<element> elements = element_map_rectangle(nbox-1, nbox-1);

	unsigned num_gp = 2*2*elements.size();
	particle_tl_weak **gauss_points = new particle_tl_weak*[num_gp];
	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i] = new particle_tl_weak(i);
	}

	gauss_int gi(elements, 2);
	gi.update_gauss_points(particles, gauss_points);

	boundary_conditions<particle_tl_weak> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_tl);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_tl);

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->quad_weight = 1.;	//uniform unit weights for fem nodes
		particles[i]->m = dx*dx*rho0;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i]->rho = rho0;
	}

	physical_constants physical_constants(nu, E, rho0);
	simulation_data sim_data(physical_constants, correction_constants());

	body<particle_tl_weak> b(particles, n, sim_data, dt, bc, gauss_points, num_gp);

	b.add_action(&derive_displacement_weak);
	b.add_action(&derive_velocity_tl_weak);
	b.add_action(&contmech_velocity_gradient_tl_gp);
	b.add_action(&contmech_continuity_gp);
	b.add_action(&material_eos_gp);
	b.add_action(&material_stress_rate_jaumann_gp);
	b.add_action(&contmech_cauchy_to_nominal_gp);
	b.add_action(&derive_stress_weak_tl);
	b.add_action(&contmech_momentum_tl_weak);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_tl_weak>>({b});
}

std::vector<body<particle_ul_weak>> tensile_ul_weak(unsigned int nbox) {
	// material constants
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double L = 1.;

	double dx = L/(nbox-1);
	double hdx = 1.;

	double dt = 1e-6;

	unsigned int np = nbox*nbox;
	unsigned int ns = (nbox-1)*(nbox-1);

	particle_ul_weak **particles = new particle_ul_weak*[np];
	particle_ul_weak **stress_points = new particle_ul_weak*[ns];

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particle_ul_weak *p = new particle_ul_weak(part_iter);

			p->x = px;
			p->y = py;

			particles[part_iter] = p;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}
			part_iter++;
		}
	}

	unsigned int sp_iter = 0;
	for (unsigned int i = 0; i < nbox-1; i++) {
		for (unsigned int j = 0; j < nbox-1; j++) {
			double px = i*dx + dx/2; double py = j*dx + dx/2;
			particle_ul_weak *p = new particle_ul_weak(sp_iter);

			p->x = px;
			p->y = py;

			stress_points[sp_iter] = p;

			sp_iter++;
		}
	}

	boundary_conditions<particle_ul_weak> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_ul);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_ul);

	for (unsigned int i = 0; i < np; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	for (unsigned int i = 0; i < ns; i++) {
		stress_points[i]->rho = rho0;
		stress_points[i]->h   = hdx*dx;
		stress_points[i]->m   = dx*dx*rho0;
		stress_points[i]->quad_weight = stress_points[i]->m/stress_points[i]->rho;
	}

	find_neighbors(particles, np, 1.5*dx, distance_manhattan);
	find_neighbors(particles, np, stress_points, ns, 1.5*dx, distance_manhattan);

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(), 0, true));

	body<particle_ul_weak> b(particles, np, sim_data, dt, bc, stress_points, ns);

	b.add_action(&riemann_update_integration_weights);
	b.add_action(&riemann_update_integration_weights_gp);

	b.add_action(&material_eos_gp);
	b.add_action(&derive_stress_weak_ul);
	b.add_action(&derive_velocity_weak);
	b.add_action(&derive_velocity);
	b.add_action(&material_stress_rate_jaumann_gp);
	b.add_action(&contmech_momentum_ul);
	b.add_action(&contmech_continuity);
	b.add_action(&contmech_continuity_gp);
	b.add_action(&contmech_advection);
	b.add_action(&interpolation_interpolate_velocity);
	b.add_action(&contmech_advection_gp);

	return std::vector<body<particle_ul_weak>>({b});
}

std::vector<body<particle_tl_weak>> tensile_tl_weak(unsigned int nbox) {
	// material constants
	double E = 1e7;
	double nu = 0.4;
	double rho0 = 1.;

	//problem dimensions
	double L = 1.;

	double dx = L/(nbox-1);
	double hdx = 1.3;

	double dt = 1e-7;

	particle_tl_weak **particles = new particle_tl_weak*[nbox*nbox];
	unsigned n = nbox*nbox;

	std::vector<unsigned int> left_bndry;
	std::vector<unsigned int> right_bndry;

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = i*dx; double py = j*dx;
			particle_tl_weak *p = new particle_tl_weak(part_iter);

			p->x = px;
			p->y = py;

			p->X = px;
			p->Y = py;

			particles[part_iter] = p;

			if (i==0) {
				left_bndry.push_back(part_iter);
			}

			if (i==nbox-1) {
				right_bndry.push_back(part_iter);
			}
			part_iter++;
		}
	}

	std::vector<element> elements = element_map_rectangle(nbox-1, nbox-1);

	unsigned num_gp = 2*2*elements.size();
	particle_tl_weak **gauss_points = new particle_tl_weak*[num_gp];
	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i] = new particle_tl_weak(i);
	}

	gauss_int gi(elements, 2);
	gi.update_gauss_points(particles, gauss_points);

	boundary_conditions<particle_tl_weak> bc;
	bc.add_boundary_condition(left_bndry,  &left_bndry_fun_tl);
	bc.add_boundary_condition(right_bndry, &right_bndry_fun_tl);

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points[i]->h   = hdx*dx;
		gauss_points[i]->rho = rho0;
	}

	find_neighbors(particles,  n, hdx*dx, distance_euclidian);
	find_neighbors(particles, n, gauss_points, num_gp, hdx*dx, distance_euclidian);

	precomp_rkpm<particle_tl_weak>(particles, gauss_points, num_gp);
	precomp_rkpm<particle_tl_weak>(particles, n);

	physical_constants physical_constants(nu, E, rho0);
	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(), 0, true));

	body<particle_tl_weak> b(particles, n, sim_data, dt, bc, gauss_points, num_gp);

	b.add_action(&derive_displacement_weak);
	b.add_action(&derive_velocity_tl_weak);
	b.add_action(&contmech_velocity_gradient_tl_gp);
	b.add_action(&contmech_continuity_gp);
	b.add_action(&material_eos_gp);
	b.add_action(&material_stress_rate_jaumann_gp);
	b.add_action(&contmech_cauchy_to_nominal_gp);
	b.add_action(&derive_stress_weak_tl);
	b.add_action(&contmech_momentum_tl_weak);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_tl_weak>>({b});
}
