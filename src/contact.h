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

//the contact algorithm in two versions:
//	contact_apply_body_to_boundary <= tests all particles in slave body for contact
//	contact_apply_boundary_to_boundary <= tests only boundary slave particles for contact (not used)
//
//sketch of contact algorithm
//	boundary of master body is sorted into a kd-tree
//	for all slave nodes:
//		extract closest point on master boundary using the kd-tree (log N query)
//		apply boundary force
//
//boundary force computation is wide spread, but for example used in:
//	Nianfei, Gan, Li Guangyao, and Long Shuyao. "3D adaptive RKPM method for contact problems with elastic–plastic dynamic large deformation." Engineering analysis with boundary elements 33.10 (2009): 1211-1222.

#ifndef CONTACT_H_
#define CONTACT_H_

#include "body.h"

template<typename T, typename Q>
void contact_apply_body_to_boundary(body<T> &master, body<Q> &slave) {
	if (!(slave.has_boundary() && master.has_boundary())) return;

	const boundary<T> *master_boundary = master.get_boundary();

	simulation_time *time = &simulation_time::getInstance();
	double dt = time->get_dt();
	double dt2 = dt*dt;

	kd_tree tree(master_boundary->points(), master_boundary->num_bound());

	const double alpha = 0.1;

	for (unsigned int i = 0; i < slave.get_num_part(); i++) {
		slave.get_particles()[i]->fcx = 0.;
		slave.get_particles()[i]->fcy = 0.;
	}

	for (unsigned int i = 0; i < slave.get_num_part(); i++) {
		kd_point query_point;
		query_point.x = slave.get_particles()[i]->x;
		query_point.y = slave.get_particles()[i]->y;

		kd_point closest_point = tree.find_closest_point(query_point.x, query_point.y);

		double dx = query_point.x - closest_point.x;
		double dy = query_point.y - closest_point.y;

		double pdotn = dx*closest_point.nx + dy*closest_point.ny;
		bool inside = pdotn <= 0.;
		if (!inside) continue;

		double gN_nx = closest_point.nx*pdotn;
		double gN_ny = closest_point.ny*pdotn;

		double ms = slave.get_particles()[i]->m;

		double fcx = ms*gN_nx/dt2*alpha;
		double fcy = ms*gN_ny/dt2*alpha;

		slave.get_particles()[i]->fcx = -fcx;
		slave.get_particles()[i]->fcy = -fcy;
	}
}

template<typename T, typename Q>
void contact_apply_boundary_to_boundary(body<T> &master, body<Q> &slave) {
	if (!(slave.has_boundary() && master.has_boundary())) return;

	const boundary<T> *master_boundary = master.get_boundary();
	const boundary<Q> *slave_boundary  = slave.get_boundary();

	simulation_time *time = &simulation_time::getInstance();
	double dt = time->get_dt();
	double dt2 = dt*dt;

	kd_tree tree(master_boundary->points(), master_boundary->num_bound());

	const double alpha = 0.1;

	for (unsigned int i = 0; i < slave.get_num_part(); i++) {
		slave.get_particles()[i]->fcx = 0.;
		slave.get_particles()[i]->fcy = 0.;
	}

	for (unsigned int i = 0; i < slave_boundary->num_bound(); i++) {
		kd_point query_point = slave_boundary->points()[i];
		kd_point closest_point = tree.find_closest_point(query_point.x, query_point.y);

		double dx = query_point.x - closest_point.x;
		double dy = query_point.y - closest_point.y;

		double pdotn = dx*closest_point.nx + dy*closest_point.ny;
		bool inside = pdotn <= 0.;
		std::vector<body<particle_ul_weak>> bodies(2);
		if (!inside) continue;

		double gN_nx = closest_point.nx*pdotn;
		double gN_ny = closest_point.ny*pdotn;

		double ms = slave.get_particles()[slave_boundary->indexes()[i]]->m;

		double fcx = ms*gN_nx/dt2*alpha;
		double fcy = ms*gN_ny/dt2*alpha;

		slave.get_particles()[slave_boundary->indexes()[i]]->fcx = -fcx;
		slave.get_particles()[slave_boundary->indexes()[i]]->fcy = -fcy;
	}
}

#endif /* CONTACT_H_ */
