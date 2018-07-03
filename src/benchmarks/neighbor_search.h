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

//some methods to identify neighbors in n2 time
//useful for total lagrangian methods where neighbors need to be computed only once

#ifndef NEIGHBOR_SEARCH_H_
#define NEIGHBOR_SEARCH_H_

#include "../particle.h"

#include <assert.h>

enum distance_type {
	distance_manhattan, distance_euclidian
};

static bool distance_r2(double xi, double yi, double xj, double yj, double h) {
	double xij = xi-xj;
	double yij = yi-yj;

	double r2 = xij*xij + yij*yij;

	double radius2 = h*h*2*2;

	return r2 <= radius2 + 1e-6;
}

static bool distance_mthn(double xi, double yi, double xj, double yj, double h) {
	double xij = xi-xj;
	double yij = yi-yj;

	return fabs(xij) < 2.*h + 1e-6 and fabs(yij) < 2.*h + 1e-6 ;
}

template <typename T>
typename std::enable_if<!std::is_base_of<particle_ul_weak, T>::value>::type
find_neighbors(T** particles, unsigned int num_part, T**gauss_points, unsigned int num_gp, double h, distance_type type) {
	for (unsigned int i = 0; i < num_gp; i++) {
		T *p1 = gauss_points[i];
		unsigned int nbh_iter = 0;

		double hi = h;

		double xi = p1->x;
		double yi = p1->y;

		for (unsigned int j = 0; j < num_part; j++) {
			T  *p2 = particles[j];

			double xj = p2->x;
			double yj = p2->y;

			bool in = (type == distance_euclidian) ? distance_r2(xi, yi, xj, yj, hi) : distance_mthn(xi, yi, xj, yj, hi);

			if (in) {
				gauss_points[i]->nbh[nbh_iter] = j;
				nbh_iter++;
			}
		}

		assert(nbh_iter < MAX_NBH);
		assert(nbh_iter > 0);

		gauss_points[i]->num_nbh = nbh_iter;
	}
}

template <typename T>
typename std::enable_if<std::is_base_of<particle_ul_weak, T>::value>::type
find_neighbors(T** particles, unsigned int num_part, T**gauss_points, unsigned int num_gp, double h, distance_type type) {
	for (unsigned int i = 0; i < num_gp; i++) {
		T *p1 = gauss_points[i];
		unsigned int nbh_iter = 0;

		double hi = h;

		double xi = p1->x;
		double yi = p1->y;

		for (unsigned int j = 0; j < num_part; j++) {
			T  *p2 = particles[j];

			double xj = p2->x;
			double yj = p2->y;

			bool in = (type == distance_euclidian) ? distance_r2(xi, yi, xj, yj, hi) : distance_mthn(xi, yi, xj, yj, hi);

			if (in) {
				gauss_points[i]->cand_nbh[nbh_iter] = j;
				nbh_iter++;
			}
		}

		assert(nbh_iter < MAX_NBH);
		assert(nbh_iter > 0);

		gauss_points[i]->num_cand_nbh = nbh_iter;
	}
}

template <typename T>
typename std::enable_if<!std::is_base_of<particle_ul_weak, T>::value>::type
find_neighbors(T** particles, unsigned int num_part, double h, distance_type type) {
	for (unsigned int i = 0; i < num_part; i++) {
		T *p1 = particles[i];
		unsigned int nbh_iter = 0;

		double hi = h;

		double xi = p1->x;
		double yi = p1->y;

		for (unsigned int j = 0; j < num_part; j++) {
			T *p2 = particles[j];

			double xj = p2->x;
			double yj = p2->y;

			bool in = (type == distance_euclidian) ? distance_r2(xi, yi, xj, yj, hi) : distance_mthn(xi, yi, xj, yj, hi);

			if (in) {
				particles[i]->nbh[nbh_iter] = j;
				nbh_iter++;
			}
		}

		assert(nbh_iter < MAX_NBH);
		assert(nbh_iter > 0);

		particles[i]->num_nbh = nbh_iter;
	}
}

template <typename T>
typename std::enable_if<std::is_base_of<particle_ul_weak, T>::value>::type
find_neighbors(T** particles, unsigned int num_part, double h, distance_type type) {
	for (unsigned int i = 0; i < num_part; i++) {
		T *p1 = particles[i];
		unsigned int nbh_iter = 0;

		double hi = h;

		double xi = p1->x;
		double yi = p1->y;

		for (unsigned int j = 0; j < num_part; j++) {
			T *p2 = particles[j];

			double xj = p2->x;
			double yj = p2->y;

			bool in = (type == distance_euclidian) ? distance_r2(xi, yi, xj, yj, hi) : distance_mthn(xi, yi, xj, yj, hi);

			if (in) {
				particles[i]->cand_nbh[nbh_iter] = j;
				nbh_iter++;
			}
		}

		assert(nbh_iter < MAX_NBH);
		assert(nbh_iter > 0);

		particles[i]->num_cand_nbh = nbh_iter;
	}
}

#endif /* NEIGHBOR_SEARCH_H_ */
