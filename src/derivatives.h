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

//methods to compute spatial derivatives different reference frames (tl, ul) and formulations (weak, strong)

#ifndef DERIVATIVES_H_
#define DERIVATIVES_H_

#include <math.h>
#include <stdio.h>
#include <glm/glm.hpp>

#include "kernel.h"
#include "particle.h"

template <typename T> class body;

template<typename T>
void derive_velocity(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double vxi = pi->vx;
		double vyi = pi->vy;

		double vx_x = 0.;
		double vx_y = 0.;
		double vy_x = 0.;
		double vy_y = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj = pj->vx;
			double vyj = pj->vy;

			double quad_weight = pj->quad_weight;

			vx_x += (vxj-vxi)*w.w_x*quad_weight;
			vx_y += (vxj-vxi)*w.w_y*quad_weight;
			vy_x += (vyj-vyi)*w.w_x*quad_weight;
			vy_y += (vyj-vyi)*w.w_y*quad_weight;
		}

		pi->vx_x = vx_x;
		pi->vx_y = vx_y;
		pi->vy_x = vy_x;
		pi->vy_y = vy_y;
	}
}

template<typename T>
void derive_velocity_weak(body<T> &b) {
	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double vx_x = 0.;
		double vx_y = 0.;
		double vy_x = 0.;
		double vy_y = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj = pj->vx;
			double vyj = pj->vy;

			double quad_weight = pj->quad_weight;

			vx_x += vxj*w.w_x*quad_weight;
			vx_y += vxj*w.w_y*quad_weight;
			vy_x += vyj*w.w_x*quad_weight;
			vy_y += vyj*w.w_y*quad_weight;
		}

		pi->vx_x = vx_x;
		pi->vx_y = vx_y;
		pi->vy_x = vy_x;
		pi->vy_y = vy_y;
	}
}

template<typename T>
void derive_displacement(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Hxx = 0.;
		double Hxy = 0.;
		double Hyx = 0.;
		double Hyy = 0.;

		double xi = pi->x;
		double yi = pi->y;

		double Xi = pi->X;
		double Yi = pi->Y;

		double uxi = xi-Xi;
		double uyi = yi-Yi;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double xj = pj->x;
			double yj = pj->y;

			double Xj = pj->X;
			double Yj = pj->Y;

			double uxj = xj-Xj;
			double uyj = yj-Yj;

			double quad_weight = pj->quad_weight;

			Hxx += (uxj-uxi)*w.w_x*quad_weight;
			Hxy += (uxj-uxi)*w.w_y*quad_weight;
			Hyx += (uyj-uyi)*w.w_x*quad_weight;
			Hyy += (uyj-uyi)*w.w_y*quad_weight;
		}

		pi->Hxx = Hxx;
		pi->Hxy = Hxy;
		pi->Hyx = Hyx;
		pi->Hyy = Hyy;
	}
}


template<typename T>
void derive_displacement_weak(body<T> &b) {
	static_assert(std::is_base_of<particle_tl_weak, T>::value, "T ist not derived from particle_tl");

	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double Hxx = 0.;
		double Hxy = 0.;
		double Hyx = 0.;
		double Hyy = 0.;

		double xi = pi->x;
		double yi = pi->y;

		double Xi = pi->X;
		double Yi = pi->Y;

		double uxi = xi-Xi;
		double uyi = yi-Yi;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double xj = pj->x;
			double yj = pj->y;

			double Xj = pj->X;
			double Yj = pj->Y;

			double uxj = xj-Xj;
			double uyj = yj-Yj;

			double quad_weight = pj->quad_weight;

			Hxx += (uxj-uxi)*w.w_x*quad_weight;
			Hxy += (uxj-uxi)*w.w_y*quad_weight;
			Hyx += (uyj-uyi)*w.w_x*quad_weight;
			Hyy += (uyj-uyi)*w.w_y*quad_weight;
		}

		pi->Hxx = Hxx;
		pi->Hxy = Hxy;
		pi->Hyx = Hyx;
		pi->Hyy = Hyy;
	}
}

template<typename T>
void derive_displacement_corotational(body<T> &b) {
	static_assert(std::is_base_of<particle_corotational, T>::value, "T ist not derived from particle_corotational");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Hxx = 0.;
		double Hxy = 0.;
		double Hyx = 0.;
		double Hyy = 0.;

		double xi = pi->x;
		double yi = pi->y;

		double Xi = pi->X;
		double Yi = pi->Y;

		double Rixx_inv = pi->Rxx;
		double Rixy_inv = pi->Ryx;
		double Riyx_inv = pi->Rxy;
		double Riyy_inv = pi->Ryy;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double xj = pj->x;
			double yj = pj->y;

			double Xj = pj->X;
			double Yj = pj->Y;

			double quad_weight = pj->quad_weight;

			double xji = xj-xi;
			double yji = yj-yi;

			double Xji = Xj-Xi;
			double Yji = Yj-Yi;

			double uxji = Rixx_inv*xji + Rixy_inv*yji - Xji;
			double uyji = Riyx_inv*xji + Riyy_inv*yji - Yji;

			Hxx += uxji*w.w_x*quad_weight;
			Hxy += uxji*w.w_y*quad_weight;
			Hyx += uyji*w.w_x*quad_weight;
			Hyy += uyji*w.w_y*quad_weight;
		}

		pi->Hxx = Hxx;
		pi->Hxy = Hxy;
		pi->Hyx = Hyx;
		pi->Hyy = Hyy;
	}
}

template<typename T>
void derive_stress_weak_tl(body<T> &b) {
	static_assert(std::is_base_of<particle_tl_weak, T>::value, "T ist not derived from particle_weak_tl");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		b.get_cur_particles()[i]->Pxx_x = 0.;
		b.get_cur_particles()[i]->Pyx_y = 0.;
		b.get_cur_particles()[i]->Pxy_x = 0.;
		b.get_cur_particles()[i]->Pyy_y = 0.;
	}

	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double Pxxi = pi->Pxx;
		double Pxyi = pi->Pxy;
		double Pyxi = pi->Pyx;
		double Pyyi = pi->Pyy;

		double quad_weight = pi->quad_weight;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vol_node = pj->quad_weight;

			pj->Pxx_x -= Pxxi*w.w_x*quad_weight*vol_node;
			pj->Pyx_y -= Pyxi*w.w_y*quad_weight*vol_node;
			pj->Pxy_x -= Pxyi*w.w_x*quad_weight*vol_node;
			pj->Pyy_y -= Pyyi*w.w_y*quad_weight*vol_node;
		}
	}
}

template<typename T>
void derive_stress_weak_ul(body<T> &b) {
	static_assert(std::is_base_of<particle_ul_weak, T>::value, "T ist not derived from particle_weak_ul");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		b.get_cur_particles()[i]->Sxx_x = 0.;
		b.get_cur_particles()[i]->Sxy_y = 0.;
		b.get_cur_particles()[i]->Sxy_x = 0.;
		b.get_cur_particles()[i]->Syy_y = 0.;
	}

	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double Sxxi = pi->Sxx - pi->p;
		double Sxyi = pi->Sxy;
		double Syyi = pi->Syy - pi->p;

		double quad_weight = pi->quad_weight;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vol_node = pj->quad_weight;

			pj->Sxx_x -= Sxxi*w.w_x*quad_weight*vol_node;
			pj->Sxy_y -= Sxyi*w.w_y*quad_weight*vol_node;
			pj->Sxy_x -= Sxyi*w.w_x*quad_weight*vol_node;
			pj->Syy_y -= Syyi*w.w_y*quad_weight*vol_node;
		}
	}
}

template<typename T>
void derive_stress_tl(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	double rho0 = b.get_sim_data().get_physical_constants().rho0();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Pxxi = pi->Pxx;
		double Pxyi = pi->Pxy;
		double Pyxi = pi->Pyx;
		double Pyyi = pi->Pyy;

		double rhoi21 = 1./(rho0*rho0);

		double Pxx_x = 0.;
		double Pyx_y = 0.;
		double Pxy_x = 0.;
		double Pyy_y = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double mj    = pj->m;
			double rhoj21 = 1./(rho0*rho0);

			double Pxxj = pj->Pxx;
			double Pxyj = pj->Pxy;
			double Pyxj = pj->Pyx;
			double Pyyj = pj->Pyy;

			Pxx_x += mj*(Pxxi*rhoi21 + Pxxj*rhoj21)*w.w_x*rho0;
			Pyx_y += mj*(Pyxi*rhoi21 + Pyxj*rhoj21)*w.w_y*rho0;
			Pxy_x += mj*(Pxyi*rhoi21 + Pxyj*rhoj21)*w.w_x*rho0;
			Pyy_y += mj*(Pyyi*rhoi21 + Pyyj*rhoj21)*w.w_y*rho0;
		}

		pi->Pxx_x = Pxx_x;
		pi->Pyx_y = Pyx_y;
		pi->Pxy_x = Pxy_x;
		pi->Pyy_y = Pyy_y;
	}
}

template<typename T>
void derive_stress_monaghan(body<T> &b) {
	static_assert(std::is_base_of<particle_monaghan, T>::value, "T ist not derived from particle_monaghan");

	double wdeltap  = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_wdeltap();
	double corr_exp = b.get_sim_data().get_correction_constants().get_monaghan_const().mghn_corr_exp();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Sxxi = pi->Sxx - pi->p;
		double Sxyi = pi->Sxy;
		double Syyi = pi->Syy - pi->p;

		double Rxxi = pi->Rxx;
		double Rxyi = pi->Rxy;
		double Ryyi = pi->Ryy;

		double rhoi = pi->rho;
		double rhoi21 = 1./(rhoi*rhoi);

		double Sxx_x = 0.;
		double Sxy_y = 0.;
		double Sxy_x = 0.;
		double Syy_y = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double Sxxj = pj->Sxx - pj->p;
			double Sxyj = pj->Sxy;
			double Syyj = pj->Syy - pj->p;

			double Rxxj = pj->Rxx;
			double Rxyj = pj->Rxy;
			double Ryyj = pj->Ryy;

			double mj     = pj->m;
			double rhoj21 = 1./(pj->rho*pj->rho);

			double Rxx = 0;
			double Rxy = 0.;
			double Ryy = 0.;

			if (wdeltap > 0) {
				double fab = w.w/wdeltap;
				fab = pow(fab,corr_exp);
				Rxx = fab*(Rxxi + Rxxj);
				Rxy = fab*(Rxyi + Rxyj);
				Ryy = fab*(Ryyi + Ryyj);
			}

			Sxx_x += mj*(Sxxi*rhoi21 + Sxxj*rhoj21 + Rxx)*w.w_x;
			Sxy_y += mj*(Sxyi*rhoi21 + Sxyj*rhoj21 + Rxy)*w.w_y;
			Sxy_x += mj*(Sxyi*rhoi21 + Sxyj*rhoj21 + Rxy)*w.w_x;
			Syy_y += mj*(Syyi*rhoi21 + Syyj*rhoj21 + Ryy)*w.w_y;
		}

		pi->Sxx_x = Sxx_x*rhoi;
		pi->Sxy_y = Sxy_y*rhoi;
		pi->Sxy_x = Sxy_x*rhoi;
		pi->Syy_y = Syy_y*rhoi;
	}
}

template<typename T>
void derive_velocity_tl_weak(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	for (unsigned int i = 0; i < b.get_num_quad_points(); i++) {
		T *pi = b.get_cur_quad_points()[i];

		double Fdotxx = 0.;
		double Fdotxy = 0.;
		double Fdotyx = 0.;
		double Fdotyy = 0.;

		double vxi = pi->vx;
		double vyi = pi->vy;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj = pj->vx;
			double vyj = pj->vy;

			double quad_weight = pj->quad_weight;

			Fdotxx += (vxj-vxi)*w.w_x*quad_weight;
			Fdotxy += (vxj-vxi)*w.w_y*quad_weight;
			Fdotyx += (vyj-vyi)*w.w_x*quad_weight;
			Fdotyy += (vyj-vyi)*w.w_y*quad_weight;
		}

		pi->Fdotxx = Fdotxx;
		pi->Fdotxy = Fdotxy;
		pi->Fdotyx = Fdotyx;
		pi->Fdotyy = Fdotyy;
	}
}

template<typename T>
void derive_velocity_tl(body<T> &b) {
	static_assert(std::is_base_of<particle_tl, T>::value, "T ist not derived from particle_tl");

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Fdotxx = 0.;
		double Fdotxy = 0.;
		double Fdotyx = 0.;
		double Fdotyy = 0.;

		double vxi = pi->vx;
		double vyi = pi->vy;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double vxj = pj->vx;
			double vyj = pj->vy;

			double quad_weight = pj->quad_weight;

			Fdotxx += (vxj-vxi)*w.w_x*quad_weight;
			Fdotxy += (vxj-vxi)*w.w_y*quad_weight;
			Fdotyx += (vyj-vyi)*w.w_x*quad_weight;
			Fdotyy += (vyj-vyi)*w.w_y*quad_weight;
		}

		pi->Fdotxx = Fdotxx;
		pi->Fdotxy = Fdotxy;
		pi->Fdotyx = Fdotyx;
		pi->Fdotyy = Fdotyy;
	}
}

template<typename T>
void derive_stress_corotational(body<T> &b) {
	static_assert(std::is_base_of<particle_corotational, T>::value, "T ist not derived from particle_corotational");

	double rho0 = b.get_sim_data().get_physical_constants().rho0();

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];

		double Pxxi = pi->Pxx;
		double Pxyi = pi->Pxy;
		double Pyxi = pi->Pyx;
		double Pyyi = pi->Pyy;

		double Rxxi = pi->Rxx;
		double Rxyi = pi->Rxy;
		double Ryxi = pi->Ryx;
		double Ryyi = pi->Ryy;

		double rhoi21 = 1./(rho0*rho0);

		double vx_t = 0.;
		double vy_t = 0.;

		for (unsigned int j = 0; j < pi->num_nbh; j++) {
			unsigned int jdx = pi->nbh[j];
			kernel_result w = pi->w[j];
			T *pj = b.get_cur_particles()[jdx];

			double Pxxj = pj->Pxx;
			double Pxyj = pj->Pxy;
			double Pyxj = pj->Pyx;
			double Pyyj = pj->Pyy;

			double Rxxj = pj->Rxx;
			double Rxyj = pj->Rxy;
			double Ryxj = pj->Ryx;
			double Ryyj = pj->Ryy;

			double mj     = pj->m;
			double rhoj21 = 1./(rho0*rho0);

			double aijx = mj*((Pxxi*rhoi21 + Pxxj*rhoj21)*w.w_x + (Pxyi*rhoi21 + Pxyj*rhoj21)*w.w_y);		//TODO: times rhoi?
			double aijy = mj*((Pyxi*rhoi21 + Pyxj*rhoj21)*w.w_x + (Pyyi*rhoi21 + Pyyj*rhoj21)*w.w_y);

			double ajix = -aijx;
			double ajiy = -aijy;

			double aijxprime = aijx*Rxxi + aijy*Rxyi;
			double aijyprime = aijx*Ryxi + aijy*Ryyi;

			double ajixprime = ajix*Rxxj + ajiy*Rxyj;
			double ajiyprime = ajix*Ryxj + ajiy*Ryyj;

			vx_t += 0.5*(aijxprime-ajixprime);
			vy_t += 0.5*(aijyprime-ajiyprime);
		}

		pi->vx_t = vx_t;
		pi->vy_t = vy_t;
	}

	for (unsigned int i = 0; i < b.get_num_part(); i++) {
		T *pi = b.get_cur_particles()[i];
		pi->vx_t += pi->fcx/pi->m;
		pi->vy_t += pi->fcy/pi->m;
	}
}

#endif /* DERIVATIVES_H_ */
