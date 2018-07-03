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

#include "rings.h"

std::vector<body<particle_monaghan>> rings_monaghan(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions (monaghan & gray)
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-6;

	particle_monaghan **particles = new particle_monaghan*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles[part_iter] = new particle_monaghan(part_iter);
				particles[part_iter]->x = px-spacing;
				particles[part_iter]->y = py;
				part_iter++;
				particles[part_iter] = new particle_monaghan(part_iter);
				particles[part_iter]->x = px+spacing;
				particles[part_iter]->y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
		particles[i]->vx = (particles[i]->x < 0.) ? 180 : -180;
	}

	//correction constants (monaghan & gray)
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

	body<particle_monaghan> b(particles, n, sim_data, dt);

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

std::vector<body<particle_ul>> rings_godunov(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-7;

	particle_ul **particles = new particle_ul*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles[part_iter] = new particle_ul(part_iter);
				particles[part_iter]->x = px-spacing;
				particles[part_iter]->y = py;
				part_iter++;
				particles[part_iter] = new particle_ul(part_iter);
				particles[part_iter]->x = px+spacing;
				particles[part_iter]->y = py;
				part_iter++;
			}
		}
	}

	unsigned int n = part_iter;

	for (unsigned int i = 0; i < n; i++) {
		particles[i]->rho = rho0;
		particles[i]->h = hdx*dx;
		particles[i]->m = dx*dx*rho0;
		particles[i]->quad_weight = particles[i]->m/particles[i]->rho;
		particles[i]->vx = (particles[i]->x < 0.) ? 180 : -180;
	}

	simulation_data sim_data(physical_constants, correction_constants());

	body<particle_ul> b(particles, n, sim_data, dt);

	b.add_action(&riemann_update_integration_weights);
	b.add_action(&material_eos);
	b.add_action(&derive_velocity);
	b.add_action(&material_stress_rate_jaumann);
	b.add_action(&godunov_continuity);
	b.add_action(&godunov_momentum);
	b.add_action(&contmech_advection);

	return std::vector<body<particle_ul>>({b});
}

std::vector<body<particle_ul>> rings_godunov_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-7;

	particle_ul **particles_left  = new particle_ul*[nbox*nbox];
	particle_ul **particles_right = new particle_ul*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles_left[part_iter] = new particle_ul(part_iter);
				particles_left[part_iter]->x = px-spacing;
				particles_left[part_iter]->y = py;

				particles_right[part_iter] = new particle_ul(part_iter);
				particles_right[part_iter]->x = px+spacing;
				particles_right[part_iter]->y = py;

				part_iter++;
			}
		}
	}

	unsigned int n_per_ring = part_iter;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx = -180.;
	}

	simulation_data sim_data(physical_constants, correction_constants());

	std::vector<body<particle_ul>> bodies(2);
	bodies[0] = body<particle_ul>(particles_left,  n_per_ring, sim_data, dt);
	bodies[1] = body<particle_ul>(particles_right, n_per_ring, sim_data, dt);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_ul> *bnd = new boundary<particle_ul>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&riemann_update_integration_weights);
		b->add_action(&material_eos);
		b->add_action(&derive_velocity);
		b->add_action(&material_stress_rate_jaumann);
		b->add_action(&godunov_continuity);
		b->add_action(&godunov_momentum);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_corotational>> rings_corot_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-7;

	particle_corotational **particles_left  = new particle_corotational*[nbox*nbox];
	particle_corotational **particles_right = new particle_corotational*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles_left[part_iter] = new particle_corotational(part_iter);
				particles_left[part_iter]->x = px - spacing;
				particles_left[part_iter]->y = py;

				particles_left[part_iter]->X = px - spacing;
				particles_left[part_iter]->Y = py;

				particles_right[part_iter] = new particle_corotational(part_iter);
				particles_right[part_iter]->x = px + spacing;
				particles_right[part_iter]->y = py;

				particles_right[part_iter]->X = px + spacing;
				particles_right[part_iter]->Y = py;

				part_iter++;
			}
		}
	}

	unsigned int n_per_ring = part_iter;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx =  -180;
	}

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(1., 1., 0.1), 0., false));

	find_neighbors(particles_left, n_per_ring, hdx*dx, distance_euclidian);
	precomp_sph<particle_corotational>(particles_left,  n_per_ring);

	find_neighbors(particles_right, n_per_ring, hdx*dx, distance_euclidian);;
	precomp_sph<particle_corotational>(particles_right, n_per_ring);

	std::vector<body<particle_corotational>> bodies(2);
	bodies[0] = body<particle_corotational>(particles_left,  n_per_ring, sim_data, dt);
	bodies[1] = body<particle_corotational>(particles_right, n_per_ring, sim_data, dt);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_corotational> *bnd = new boundary<particle_corotational>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&corotational_anisotropy_tensor);
		b->add_action(&derive_displacement_corotational);
		b->add_action(&material_compute_stress_cauchy_green);
		b->add_action(&contmech_cauchy_to_nominal);
		b->add_action(&derive_stress_corotational);
		b->add_action(&correctors_mghn_artificial_viscosity);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_tl>> rings_tl_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-7;

	particle_tl **particles_left  = new particle_tl*[nbox*nbox];
	particle_tl **particles_right = new particle_tl*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles_left[part_iter] = new particle_tl(part_iter);
				particles_left[part_iter]->x = px - spacing;
				particles_left[part_iter]->y = py;

				particles_left[part_iter]->X = px - spacing;
				particles_left[part_iter]->Y = py;

				particles_right[part_iter] = new particle_tl(part_iter);
				particles_right[part_iter]->x = px + spacing;
				particles_right[part_iter]->y = py;

				particles_right[part_iter]->X = px + spacing;
				particles_right[part_iter]->Y = py;

				part_iter++;
			}
		}
	}

	unsigned int n_per_ring = part_iter;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx =  -180;
	}

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(3., 3., 0.1), 0., false));

	find_neighbors(particles_left, n_per_ring, hdx*dx, distance_euclidian);
	precomp_randles_libersky<particle_tl>(particles_left,  n_per_ring);

	find_neighbors(particles_right, n_per_ring, hdx*dx, distance_euclidian);;
	precomp_randles_libersky<particle_tl>(particles_right, n_per_ring);

	std::vector<body<particle_tl>> bodies(2);
	bodies[0] = body<particle_tl>(particles_left,  n_per_ring, sim_data, dt);
	bodies[1] = body<particle_tl>(particles_right, n_per_ring, sim_data, dt);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_tl> *bnd = new boundary<particle_tl>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&derive_displacement);
		b->add_action(&material_compute_stress_green_stvenant);
		b->add_action(&contmech_cauchy_to_nominal);
		b->add_action(&derive_stress_tl);
		b->add_action(&correctors_mghn_artificial_viscosity);
		b->add_action(&contmech_momentum_tl);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_potential>>    rings_potential_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.75;
	double ro = 1.00;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.3;

	double dt = 1e-7;

	particle_potential **particles_left  = new particle_potential*[nbox*nbox];
	particle_potential **particles_right = new particle_potential*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles_left[part_iter] = new particle_potential(part_iter);
				particles_left[part_iter]->x = px - spacing;
				particles_left[part_iter]->y = py;

				particles_left[part_iter]->X = px - spacing;
				particles_left[part_iter]->Y = py;

				particles_right[part_iter] = new particle_potential(part_iter);
				particles_right[part_iter]->x = px + spacing;
				particles_right[part_iter]->y = py;

				particles_right[part_iter]->X = px + spacing;
				particles_right[part_iter]->Y = py;

				part_iter++;
			}
		}
	}

	unsigned int n_per_ring = part_iter;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx =  -180;
	}

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(1., 1., 0.1), 0., false));

	find_neighbors(particles_left, n_per_ring, hdx*dx, distance_euclidian);
	precomp_randles_libersky_symmetric<particle_potential>(particles_left,  n_per_ring);

	find_neighbors(particles_right, n_per_ring, hdx*dx, distance_euclidian);
	precomp_randles_libersky_symmetric<particle_potential>(particles_right, n_per_ring);

	std::vector<body<particle_potential>> bodies(2);
	bodies[0] = body<particle_potential>(particles_left,  n_per_ring, sim_data, dt);
	bodies[1] = body<particle_potential>(particles_right, n_per_ring, sim_data, dt);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_potential> *bnd = new boundary<particle_potential>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&derive_displacement);
		b->add_action(&material_compute_stress_green_stvenant);
		b->add_action(&potential_comp_elastic_potential);
		b->add_action(&correctors_mghn_artificial_viscosity);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_monaghan>> rings_monaghan_contact(unsigned int nbox) {

	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.7;

	double dt = 1e-6;

	particle_monaghan **particles_left  = new particle_monaghan*[nbox*nbox];
	particle_monaghan **particles_right = new particle_monaghan*[nbox*nbox];

	unsigned int part_iter = 0;
	for (unsigned int i = 0; i < nbox; i++) {
		for (unsigned int j = 0; j < nbox; j++) {
			double px = -ro+i*dx; double py = -ro+j*dx;
			double dist = sqrt(px*px + py*py);
			if (dist < ro && dist >= ri) {
				particles_left[part_iter] = new particle_monaghan(part_iter);
				particles_left[part_iter]->x = px-spacing;
				particles_left[part_iter]->y = py;

				particles_right[part_iter] = new particle_monaghan(part_iter);
				particles_right[part_iter]->x = px+spacing;
				particles_right[part_iter]->y = py;

				part_iter++;
			}
		}
	}

	unsigned int n_per_ring = part_iter;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx = -180.;
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

	std::vector<body<particle_monaghan>> bodies(2);
	bodies[0] = body<particle_monaghan>(particles_left,  n_per_ring, sim_data, dt);
	bodies[1] = body<particle_monaghan>(particles_right, n_per_ring, sim_data, dt);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_monaghan> *bnd = new boundary<particle_monaghan>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&riemann_update_integration_weights);
		b->add_action(&material_eos);
		b->add_action(&correctors_mghn_artificial_stress);
		b->add_action(&derive_stress_monaghan);
		b->add_action(&derive_velocity);
		b->add_action(&correctors_mghn_artificial_viscosity);
		b->add_action(&correctors_xsph);
		b->add_action(&material_stress_rate_jaumann);
		b->add_action(&contmech_continuity);
		b->add_action(&contmech_momentum_ul);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_ul_weak>> rings_ul_fem_contact(unsigned int nbox) {

	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);

	double dt = 1e-7;

	unsigned int n_left;
	unsigned int n_right;
	std::vector<element> element_map_left;
	std::vector<element> element_map_right;
	particle_ul_weak ** particles_left  = mesh_ring<particle_ul_weak>(nbox, n_left,  element_map_left,  ri, ro);
	particle_ul_weak ** particles_right = mesh_ring<particle_ul_weak>(nbox, n_right, element_map_right, ri, ro);

	assert(n_left == n_right);
	unsigned int n_per_ring = n_left;

	unsigned int order = 2;

	gauss_int *g_left = new gauss_int(element_map_left, order);
	unsigned int num_gp_left = order*order*element_map_left.size();
	particle_ul_weak **gauss_points_left = new particle_ul_weak*[num_gp_left];
	for (unsigned int i = 0; i < num_gp_left; i++) {
		gauss_points_left[i] = new particle_ul_weak(i);
	}
	g_left->update_gauss_points(particles_left, gauss_points_left);

	gauss_int *g_right = new gauss_int(element_map_right, order);
	unsigned int num_gp_right = order*order*element_map_right.size();
	particle_ul_weak **gauss_points_right = new particle_ul_weak*[num_gp_right];
	for (unsigned int i = 0; i < num_gp_right; i++) {
		gauss_points_right[i] = new particle_ul_weak(i);
	}
	g_right->update_gauss_points(particles_right, gauss_points_right);

	assert(num_gp_left == num_gp_right);
	unsigned int num_gp = num_gp_left;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->x -= spacing;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = 1.;
		particles_left[i]->vx =  180;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->x += spacing;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = 1.;
		particles_right[i]->vx = -180.;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_left[i]->rho = rho0;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_right[i]->rho = rho0;
	}

	simulation_data sim_data(physical_constants, correction_constants());

	std::vector<body<particle_ul_weak>> bodies(2);
	bodies[0] = body<particle_ul_weak>(particles_left,  n_per_ring, sim_data, dt, boundary_conditions<particle_ul_weak>(), gauss_points_left,  num_gp, g_left);
	bodies[1] = body<particle_ul_weak>(particles_right, n_per_ring, sim_data, dt, boundary_conditions<particle_ul_weak>(), gauss_points_right, num_gp, g_right);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_ul_weak> *bnd = new boundary<particle_ul_weak>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&material_eos_gp);
		b->add_action(&derive_stress_weak_ul);
		b->add_action(&derive_velocity_weak);
		b->add_action(&material_stress_rate_jaumann_gp);
		b->add_action(&contmech_momentum_ul_weak);
		b->add_action(&contmech_continuity_gp);
		b->add_action(&contmech_advection);
	}

	std::vector<body<particle_ul_weak>> dummy;
	dummy.push_back(bodies[0]);

	return bodies;
}

std::vector<body<particle_tl_weak>> rings_tl_fem_contact(unsigned int nbox) {

	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);

	double dt = 1e-7;

	unsigned int n_left;
	unsigned int n_right;
	std::vector<element> element_map_left;
	std::vector<element> element_map_right;
	particle_tl_weak **particles_left  = mesh_ring<particle_tl_weak>(nbox, n_left,  element_map_left,  ri, ro);
	particle_tl_weak **particles_right = mesh_ring<particle_tl_weak>(nbox, n_right, element_map_right, ri, ro);

	assert(n_left == n_right);
	unsigned int n_per_ring = n_left;

	unsigned int order = 2;

	gauss_int *g_left = new gauss_int(element_map_left, order);
	unsigned int num_gp_left = order*order*element_map_left.size();
	particle_tl_weak **gauss_points_left = new particle_tl_weak*[num_gp_left];
	for (unsigned int i = 0; i < num_gp_left; i++) {
		gauss_points_left[i] = new particle_tl_weak(i);
	}
	g_left->update_gauss_points(particles_left, gauss_points_left);

	gauss_int *g_right = new gauss_int(element_map_right, order);
	unsigned int num_gp_right = order*order*element_map_right.size();
	particle_tl_weak **gauss_points_right = new particle_tl_weak*[num_gp_right];
	for (unsigned int i = 0; i < num_gp_right; i++) {
		gauss_points_right[i] = new particle_tl_weak(i);
	}
	g_right->update_gauss_points(particles_right, gauss_points_right);

	assert(num_gp_left == num_gp_right);
	unsigned int num_gp = num_gp_left;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->x -= spacing;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = 1.;
		particles_left[i]->vx = 180.;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_left[i]->rho = rho0;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->x += spacing;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = 1.;
		particles_right[i]->vx = -180.;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_right[i]->rho = rho0;
	}

	simulation_data sim_data(physical_constants, correction_constants());

	std::vector<body<particle_tl_weak>> bodies(2);
	bodies[0] = body<particle_tl_weak>(particles_left,  n_per_ring, sim_data, dt, boundary_conditions<particle_tl_weak>(), gauss_points_left,  num_gp);
	bodies[1] = body<particle_tl_weak>(particles_right, n_per_ring, sim_data, dt, boundary_conditions<particle_tl_weak>(), gauss_points_right, num_gp);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_tl_weak> *bnd = new boundary<particle_tl_weak>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&derive_displacement_weak);
		b->add_action(&derive_velocity_tl_weak);
		b->add_action(&contmech_velocity_gradient_tl_gp);
		b->add_action(&contmech_continuity_gp);
		b->add_action(&material_eos_gp);
		b->add_action(&material_stress_rate_jaumann_gp);
		b->add_action(&contmech_cauchy_to_nominal_gp);
		b->add_action(&derive_stress_weak_tl);
		b->add_action(&contmech_momentum_tl_weak);
		b->add_action(&contmech_advection);
	}

	return bodies;
}

std::vector<body<particle_ul_weak>> rings_ul_weak_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.;

	double dt = 1e-7;

	unsigned int n_left;
	unsigned int n_right;
	std::vector<element> element_map_left;
	std::vector<element> element_map_right;
	particle_ul_weak **particles_left  = mesh_ring<particle_ul_weak>(nbox, n_left,  element_map_left,  ri, ro);
	particle_ul_weak **particles_right = mesh_ring<particle_ul_weak>(nbox, n_right, element_map_right, ri, ro);

	assert(n_left == n_right);
	unsigned int n_per_ring = n_left;

	unsigned int order = 1;

	gauss_int *g_left = new gauss_int(element_map_left, order);
	unsigned int num_sp_left = order*order*element_map_left.size();
	particle_ul_weak **stress_points_left = new particle_ul_weak*[num_sp_left];
	for (unsigned int i = 0; i < num_sp_left; i++) {
		stress_points_left[i] = new particle_ul_weak(i);
	}
	g_left->update_gauss_points(particles_left, stress_points_left);

	gauss_int *g_right = new gauss_int(element_map_right, order);
	unsigned int num_sp_right = order*order*element_map_right.size();
	particle_ul_weak **stress_points_right = new particle_ul_weak*[num_sp_right];
	for (unsigned int i = 0; i < num_sp_right; i++) {
		stress_points_right[i] = new particle_ul_weak(i);
	}
	g_right->update_gauss_points(particles_right, stress_points_right);

	assert(num_sp_left == num_sp_right);
	unsigned int num_sp = num_sp_left;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->x -= spacing;
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx = 180.;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->x += spacing;
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx = -180.;
	}

	for (unsigned int i = 0; i < num_sp; i++) {
		stress_points_left[i]->x  -= spacing;
		stress_points_left[i]->rho = rho0;
		stress_points_left[i]->h   = hdx*dx;
		stress_points_left[i]->m   = dx*dx*rho0;
		stress_points_left[i]->quad_weight = stress_points_left[i]->m/stress_points_left[i]->rho;
	}

	for (unsigned int i = 0; i < num_sp; i++) {
		stress_points_right[i]->x  += spacing;
		stress_points_right[i]->rho = rho0;
		stress_points_right[i]->h   = hdx*dx;
		stress_points_right[i]->m   = dx*dx*rho0;
		stress_points_right[i]->quad_weight = stress_points_right[i]->m/stress_points_right[i]->rho;
	}

	find_neighbors(particles_left,  n_per_ring, 1.5*dx, distance_manhattan);
	find_neighbors(particles_right, n_per_ring, 1.5*dx, distance_manhattan);

	find_neighbors(particles_left,  n_per_ring, stress_points_left,  num_sp, 1.5*dx, distance_manhattan);
	find_neighbors(particles_right, n_per_ring, stress_points_right, num_sp, 1.5*dx, distance_manhattan);

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(), 0, true));

	std::vector<body<particle_ul_weak>> bodies(2);
	bodies[0] = body<particle_ul_weak>(particles_left,  n_per_ring, sim_data, dt, boundary_conditions<particle_ul_weak>(), stress_points_left,  num_sp);
	bodies[1] = body<particle_ul_weak>(particles_right, n_per_ring, sim_data, dt, boundary_conditions<particle_ul_weak>(), stress_points_right, num_sp);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_ul_weak> *bnd = new boundary<particle_ul_weak>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&riemann_update_integration_weights);
		b->add_action(&riemann_update_integration_weights_gp);

		b->add_action(&material_eos_gp);
		b->add_action(&derive_stress_weak_ul);
		b->add_action(&derive_velocity_weak);
		b->add_action(&derive_velocity);
		b->add_action(&material_stress_rate_jaumann_gp);
		b->add_action(&contmech_momentum_ul_weak);
		b->add_action(&contmech_continuity);
		b->add_action(&contmech_continuity_gp);
		b->add_action(&contmech_advection);
		b->add_action(&interpolation_interpolate_velocity);
		b->add_action(&contmech_advection_gp);
	}

	return bodies;
}

std::vector<body<particle_tl_weak>> rings_tl_weak_contact(unsigned int nbox) {
	// material constants (rubber like)
	double E    = 1e7;
	double nu   = 0.4;
	double rho0 = 1.;

	physical_constants physical_constants(nu, E, rho0);

	//problem dimensions
	double ri = 0.03;
	double ro = 0.04;
	double spacing = ro + 0.001;

	double dx = 2*ro/(nbox-1);
	double hdx = 1.;

	double dt = 1e-7;

	unsigned int n_left;
	unsigned int n_right;
	std::vector<element> element_map_left;
	std::vector<element> element_map_right;
	particle_tl_weak **particles_left  = mesh_ring<particle_tl_weak>(nbox, n_left,  element_map_left,  ri, ro);
	particle_tl_weak **particles_right = mesh_ring<particle_tl_weak>(nbox, n_right, element_map_right, ri, ro);

	assert(n_left == n_right);
	unsigned int n_per_ring = n_left;

	unsigned int order = 2;

	gauss_int *g_left = new gauss_int(element_map_left, order);
	unsigned int num_gp_left = order*order*element_map_left.size();
	particle_tl_weak **gauss_points_left = new particle_tl_weak*[num_gp_left];
	for (unsigned int i = 0; i < num_gp_left; i++) {
		gauss_points_left[i] = new particle_tl_weak(i);
	}
	g_left->update_gauss_points(particles_left, gauss_points_left);

	gauss_int *g_right = new gauss_int(element_map_right, order);
	unsigned int num_gp_right = order*order*element_map_right.size();
	particle_tl_weak **gauss_points_right = new particle_tl_weak*[num_gp_right];
	for (unsigned int i = 0; i < num_gp_right; i++) {
		gauss_points_right[i] = new particle_tl_weak(i);
	}
	g_right->update_gauss_points(particles_right, gauss_points_right);

	assert(num_gp_left == num_gp_right);
	unsigned int num_gp = num_gp_left;

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_left[i]->x -= spacing;
		particles_left[i]->rho = rho0;
		particles_left[i]->h = hdx*dx;
		particles_left[i]->m = dx*dx*rho0;
		particles_left[i]->quad_weight = particles_left[i]->m/particles_left[i]->rho;
		particles_left[i]->vx = +180.;
	}

	for (unsigned int i = 0; i < n_per_ring; i++) {
		particles_right[i]->x += spacing;
		particles_right[i]->rho = rho0;
		particles_right[i]->h = hdx*dx;
		particles_right[i]->m = dx*dx*rho0;
		particles_right[i]->quad_weight = particles_right[i]->m/particles_right[i]->rho;
		particles_right[i]->vx = -180.;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_left[i]->x  -= spacing;
		gauss_points_left[i]->h   = hdx*dx;
		gauss_points_left[i]->rho = rho0;
	}

	for (unsigned int i = 0; i < num_gp; i++) {
		gauss_points_right[i]->x  += spacing;
		gauss_points_right[i]->h   = hdx*dx;
		gauss_points_right[i]->rho = rho0;
	}

	find_neighbors(particles_left,  n_per_ring, hdx*dx, distance_euclidian);
	find_neighbors(particles_right, n_per_ring, hdx*dx, distance_euclidian);

	find_neighbors(particles_left,  n_per_ring, gauss_points_left,  num_gp, hdx*dx, distance_euclidian);
	find_neighbors(particles_right, n_per_ring, gauss_points_right, num_gp, hdx*dx, distance_euclidian);

	precomp_rkpm<particle_tl_weak>(particles_left, gauss_points_left, num_gp);
	precomp_rkpm<particle_tl_weak>(particles_left, n_per_ring);

	precomp_rkpm<particle_tl_weak>(particles_right, gauss_points_right, num_gp);
	precomp_rkpm<particle_tl_weak>(particles_right, n_per_ring);

	simulation_data sim_data(physical_constants, correction_constants(constants_monaghan(), constants_artificial_viscosity(), 0, true));

	std::vector<body<particle_tl_weak>> bodies(2);
	bodies[0] = body<particle_tl_weak>(particles_left,  n_per_ring, sim_data, dt, boundary_conditions<particle_tl_weak>(), gauss_points_left,  num_gp);
	bodies[1] = body<particle_tl_weak>(particles_right, n_per_ring, sim_data, dt, boundary_conditions<particle_tl_weak>(), gauss_points_right, num_gp);

	for (auto b = bodies.begin(); b != bodies.end(); ++b) {
		std::vector<unsigned int> boundary_idx = convex_hull(b->get_particles(), n_per_ring);
		boundary<particle_tl_weak> *bnd = new boundary<particle_tl_weak>(boundary_idx, b->get_particles());
		b->set_boundary(bnd);

		b->add_action(&derive_displacement_weak);
		b->add_action(&derive_velocity_tl_weak);
		b->add_action(&contmech_velocity_gradient_tl_gp);
		b->add_action(&contmech_continuity_gp);
		b->add_action(&material_eos_gp);
		b->add_action(&material_stress_rate_jaumann_gp);
		b->add_action(&contmech_cauchy_to_nominal_gp);
		b->add_action(&derive_stress_weak_tl);
		b->add_action(&contmech_momentum_tl_weak);
		b->add_action(&contmech_advection);
	}

	return bodies;
}
