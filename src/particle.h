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

//particle type hierarchy
//	type hierarchy to enable conditional compilation
//variables without underscore are state varibles
//variables with _t are time derivatives
//variables with _qq (q=x,y) are spatial derivatives

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "kernel.h"
#include <string.h>
#include <assert.h>

#define MAX_NBH 200

class particle {
public:
	particle();
	particle(unsigned int idx);
	virtual ~particle();

	virtual void reset();
	virtual void copy_into(particle *p) const;

	unsigned int idx = 0, hash = 0;
	double x = 0., y = 0.;
	double X = 0.; double Y = 0.;
	double m = 0., rho = 0., h = 0., quad_weight = 0.; double p = 0.;
	double vx = 0., vy = 0.;
	double Sxx = 0., Sxy = 0., Syy = 0.;

	double vx_x = 0., vx_y = 0., vy_x = 0., vy_y = 0.;

	double x_t = 0., y_t = 0.;
	double rho_t = 0.;
	double vx_t = 0., vy_t = 0.;
	double Sxx_t = 0., Sxy_t = 0., Syy_t = 0.;

	double fcx = 0.; double fcy = 0.;

	unsigned int num_nbh = 0.;
	unsigned int nbh[MAX_NBH];
	kernel_result w[MAX_NBH];


};

class particle_tl : public particle {
public:
	particle_tl();
	particle_tl(unsigned int idx);
	virtual ~particle_tl();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;

	double Hxx = 0., Hxy = 0., Hyx = 0., Hyy = 0.;
	double Fdotxx = 0.; double Fdotxy = 0.; double Fdotyx = 0.; double Fdotyy = 0.;
	double Pxx = 0., Pxy = 0., Pyx = 0., Pyy = 0.;
	double Pxx_x = 0., Pxy_x = 0., Pyx_y = 0., Pyy_y = 0.;
};

class particle_corotational : public particle_tl  {
public:
	particle_corotational();
	particle_corotational(unsigned int idx);
	virtual ~particle_corotational();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;

	double Rxx = 0., Rxy = 0., Ryx = 0., Ryy = 0.;
};

class particle_potential : public particle_tl {
public:
	particle_potential(unsigned int idx);
	particle_potential();
	virtual ~particle_potential();

	virtual void reset() override;
	void copy_into(particle * other) const override;

	kernel_result wji[MAX_NBH];
};

class particle_tl_weak : public particle_tl {
public:
	particle_tl_weak();
	particle_tl_weak(unsigned int idx);
	virtual ~particle_tl_weak();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;
};

class particle_ul : public particle {
public:
	particle_ul();
	particle_ul(unsigned int idx);
	virtual ~particle_ul();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;

	double Sxx_x = 0., Sxy_y = 0., Sxy_x = 0., Syy_y = 0.;
};

class particle_ul_weak : public particle_ul {
public:
	particle_ul_weak();
	particle_ul_weak(unsigned int idx);
	~particle_ul_weak();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;

	unsigned int num_cand_nbh = 0;
	unsigned int cand_nbh[MAX_NBH];
};

class particle_monaghan : public particle_ul {
public:
	particle_monaghan();
	particle_monaghan(unsigned int idx);
	~particle_monaghan();

	virtual void reset() override;
	virtual void copy_into(particle *p) const override;

	double Rxx = 0., Rxy = 0., Ryy = 0.;
};

#endif /* PARTICLE_H_ */
