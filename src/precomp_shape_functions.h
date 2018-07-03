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

//compute sph, rkpm, or randles libersky corrected kernels

#ifndef PRECOMPSHAPEFUNCTIONS_H_
#define PRECOMPSHAPEFUNCTIONS_H_

#include "particle.h"
#include "kernel.h"

#include <assert.h>
#include <glm/glm.hpp>

#include "gauss_integration.h"

static void solve3x3(glm::dmat3x3 M, glm::dvec3 rhs, glm::dvec3 &x) {

	double delta   =  M[0][0]*M[1][1]*M[2][2] + M[1][0]*M[2][1]*M[0][2] + M[2][0]*M[0][1]*M[1][2] -
			M[0][2]*M[1][1]*M[2][0] - M[1][2]*M[2][1]*M[0][0] - M[2][2]*M[0][1]*M[1][0];

	double delta_x =  rhs[0]*M[1][1]*M[2][2] + M[1][0]*M[2][1]*rhs[2] + M[2][0]*rhs[1]*M[1][2] -
			rhs[2]*M[1][1]*M[2][0] - M[1][2]*M[2][1]*rhs[0] - M[2][2]*rhs[1]*M[1][0];

	double delta_y =  M[0][0]*rhs[1]*M[2][2] + rhs[0]*M[2][1]*M[0][2] + M[2][0]*M[0][1]*rhs[2] -
			M[0][2]*rhs[1]*M[2][0] - rhs[2]*M[2][1]*M[0][0] - M[2][2]*M[0][1]*rhs[0];

	double delta_z =  M[0][0]*M[1][1]*rhs[2] + M[1][0]*rhs[1]*M[0][2] + rhs[0]*M[0][1]*M[1][2] -
			M[0][2]*M[1][1]*rhs[0] - M[1][2]*rhs[1]*M[0][0] - rhs[2]*M[0][1]*M[1][0];

	assert(fabs(delta) > 0.);

	x[0] = delta_x/delta;
	x[1] = delta_y/delta;
	x[2] = delta_z/delta;
}

template<typename T>
static void compute_rkpm_correctors(T * pi, T **particles, glm::dvec3 &C, glm::dvec3 &C_x, glm::dvec3 &C_y) {
	double xi = pi->x;
	double yi = pi->y;
	double hi = pi->h;

	unsigned int num_nbh = pi->num_nbh;

	glm::dmat3x3 M(0.);
	glm::dmat3x3 M_x(0.);
	glm::dmat3x3 M_y(0.);

	for (unsigned int j = 0; j < num_nbh; j++) {
		unsigned int jdx = pi->nbh[j];

		T *pj = particles[jdx];

		double mj = pj->m;
		double rhoj = pj->rho;

		double xj = pj->x;
		double yj = pj->y;

		double rx   = xi-xj;
		double ry   = yi-yj;

		double rx_x = 1.;
		double rx_y = 0.;

		double ry_x = 0.;
		double ry_y = 1.;

		kernel_result w = cubic_spline(xi,yi,xj,yj,hi);

		// moment matrix
		M[0][0] += w.w*mj/rhoj;

		M[0][1] += w.w*rx*mj/rhoj;
		M[0][2] += w.w*ry*mj/rhoj;

		M[1][0] += w.w*rx*mj/rhoj;
		M[2][0] += w.w*ry*mj/rhoj;

		M[1][1] += rx*rx*w.w*mj/rhoj;
		M[1][2] += rx*ry*w.w*mj/rhoj;
		M[2][1] += ry*rx*w.w*mj/rhoj;
		M[2][2] += ry*ry*w.w*mj/rhoj;

		//-------------------------------------------------------

		M_x[0][0] += w.w_x*mj/rhoj;

		M_x[0][1] += w.w_x*rx*mj/rhoj + w.w*rx_x*mj/rhoj;
		M_x[0][2] += w.w_x*ry*mj/rhoj + w.w*rx_y*mj/rhoj;

		M_x[1][0] += w.w_x*rx*mj/rhoj + w.w*rx_x*mj/rhoj;
		M_x[2][0] += w.w_x*ry*mj/rhoj + w.w*rx_y*mj/rhoj;

		M_x[1][1] += rx*rx*w.w_x*mj/rhoj + w.w*(rx_x*rx + rx*rx_x)*mj/rhoj;
		M_x[1][2] += rx*ry*w.w_x*mj/rhoj + w.w*(rx_x*ry + rx*rx_y)*mj/rhoj;
		M_x[2][1] += rx*ry*w.w_x*mj/rhoj + w.w*(rx_y*rx + ry*rx_x)*mj/rhoj;
		M_x[2][2] += ry*ry*w.w_x*mj/rhoj + w.w*(rx_y*ry + ry*rx_y)*mj/rhoj;

		//-------------------------------------------------------

		M_y[0][0] += w.w_y*mj/rhoj;

		M_y[0][1] += w.w_y*rx*mj/rhoj + w.w*ry_x*mj/rhoj;
		M_y[0][2] += w.w_y*ry*mj/rhoj + w.w*ry_y*mj/rhoj;

		M_y[1][0] += w.w_y*rx*mj/rhoj + w.w*ry_x*mj/rhoj;
		M_y[2][0] += w.w_y*ry*mj/rhoj + w.w*ry_y*mj/rhoj;

		M_y[1][1] += rx*rx*w.w_y*mj/rhoj + w.w*(ry_x*rx + rx*ry_x)*mj/rhoj;
		M_y[1][2] += rx*ry*w.w_y*mj/rhoj + w.w*(ry_x*ry + rx*ry_y)*mj/rhoj;
		M_y[2][1] += rx*ry*w.w_y*mj/rhoj + w.w*(ry_y*rx + ry*ry_x)*mj/rhoj;
		M_y[2][2] += ry*ry*w.w_y*mj/rhoj + w.w*(ry_y*ry + ry*ry_y)*mj/rhoj;
	}

	glm::dvec3 P(0.);
	P[0] = 1.;

	solve3x3(M,P,C);

	glm::dvec3 rhs_x = -M_x*C;
	glm::dvec3 rhs_y = -M_y*C;

	solve3x3(M,rhs_x,C_x);
	solve3x3(M,rhs_y,C_y);
}

template<typename T>
static kernel_result compute_rkpm_shape_function(T *pi, T *pj, glm::dvec3 C, glm::dvec3 C_x, glm::dvec3 C_y) {
	double xi = pi->x;
	double yi = pi->y;
	double hi = pi->h;

	double xj = pj->x;
	double yj = pj->y;

	double rx   = xi-xj;
	double ry   = yi-yj;

	double rx_x = 1.;
	double rx_y = 0.;

	double ry_x = 0.;
	double ry_y = 1.;

	kernel_result w = cubic_spline(xi,yi,xj,yj,hi);

	double Corr   = (rx*C[1]   + ry*C[2])   + C[0];
	double Corr_x = (rx*C_x[1] + ry*C_x[2]) + rx_x*C[1] + rx_y*C[2] + C_x[0];
	double Corr_y = (rx*C_y[1] + ry*C_y[2]) + ry_x*C[1] + ry_y*C[2] + C_y[0];

	kernel_result N;
	N.w   = Corr*w.w;
	N.w_x = Corr_x*w.w + Corr*w.w_x;
	N.w_y = Corr_y*w.w + Corr*w.w_y;

	return N;
}

template <typename T>
void precomp_rkpm(T **particles, unsigned int n) {
	for (unsigned int i = 0; i < n; i++) {
		T *pi = particles[i];

		glm::dvec3 C, C_x, C_y;
		compute_rkpm_correctors(pi, particles, C, C_x, C_y);

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			T *pj = particles[jdx];

			pi->w[j] = compute_rkpm_shape_function(pi, pj, C, C_x, C_y);
		}
	}
}

template <typename T>
void precomp_rkpm(T **from, T **to, unsigned int nto) {
	for (unsigned int i = 0; i < nto; i++) {
		T *pi = to[i];

		glm::dvec3 C, C_x, C_y;
		compute_rkpm_correctors(pi, from, C, C_x, C_y);

		unsigned int num_nbh = pi->num_nbh;
		for (unsigned int j = 0; j < num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			T *pj = from[jdx];

			pi->w[j] = compute_rkpm_shape_function(pi, pj, C, C_x, C_y);
		}
	}
}

template <typename T>
void precomp_sph(T **particles, unsigned int n) {
	for (unsigned int i = 0; i < n; i++) {
		particle *pi = particles[i];

		double xi = pi->x;
		double yi = pi->y;
		double hi = pi->h;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			particle *pj = particles[jdx];

			double xj = pj->x;
			double yj = pj->y;

			pi->w[j] = cubic_spline(xi, yi, xj, yj, hi);
		}
	}
}

template <typename T>
void precomp_randles_libersky(T **particles, unsigned int np) {
	for (unsigned int i = 0; i < np; i++) {
		unsigned int num_nbh = particles[i]->num_nbh;
		T *pi = particles[i];

		double xi = pi->x;
		double yi = pi->y;
		double hi = pi->h;

		double C = 0.;

		for (unsigned int j = 0; j < num_nbh; j++) {
			unsigned int jdx = particles[i]->nbh[j];
			T *pj = particles[jdx];

			double xj   = pj->x;
			double yj   = pj->y;
			double volj = pj->m/pj->rho;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

			C += w.w*volj;
		}

		glm::dmat2x2 B(0);

		for (unsigned int j = 0; j < num_nbh; j++) {
			unsigned int jdx = particles[i]->nbh[j];
			T *pj = particles[jdx];

			double xj   = pj->x;
			double yj   = pj->y;
			double volj = pj->m/pj->rho;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

			B[0][0] += (xj - xi)*w.w_x*volj*1./C;
			B[1][0] += (xj - xi)*w.w_y*volj*1./C;
			B[0][1] += (yj - yi)*w.w_x*volj*1./C;
			B[1][1] += (yj - yi)*w.w_y*volj*1./C;
		}

		glm::dmat2x2 invB = glm::inverse(B);

		for (unsigned int j = 0; j < num_nbh; j++) {
			unsigned int jdx = particles[i]->nbh[j];
			T *pj = particles[jdx];

			double xj = pj->x;
			double yj = pj->y;

			kernel_result w = cubic_spline(xi, yi, xj, yj, hi);

			pi->w[j].w = w.w*1./C;
			pi->w[j].w_x = (w.w_x*invB[0][0] + w.w_y*invB[1][0])*1./C;
			pi->w[j].w_y = (w.w_x*invB[0][1] + w.w_y*invB[1][1])*1./C;
		}
	}
}

template <typename T>
void precomp_randles_libersky_symmetric(T **particles, unsigned int np) {

	precomp_randles_libersky(particles, np);

	for (unsigned int i = 0; i < np; i++) {
		T *pi = particles[i];
		unsigned int num_nbhi = pi->num_nbh;

		for (unsigned int j = 0; j < num_nbhi; j++) {
			unsigned int jdx = pi->nbh[j];
			T *pj = particles[jdx];
			unsigned int num_nbhj = pj->num_nbh;

			for (unsigned int k = 0; k < num_nbhj; k++) {
				unsigned int kdx = pj->nbh[k];
				if (kdx == i) {
					pi->wji[j] = pj->w[k];
					break;
				}
			}
		}
	}
}

#endif /* PRECOMPSHAPEFUNCTIONS_H_ */
