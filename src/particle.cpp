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

#include "particle.h"

particle::particle() {
	memset(nbh, 0, sizeof(unsigned int)*MAX_NBH);
	memset(w,   0, sizeof(kernel_result)*MAX_NBH);
};

particle::particle(unsigned int idx) {
	this->idx = idx;

	memset(nbh, 0, sizeof(unsigned int)*MAX_NBH);
	memset(w,   0, sizeof(kernel_result)*MAX_NBH);
}

particle::~particle() {};

void particle::reset() {
	x_t  = 0.;
	y_t  = 0.;
	rho_t = 0.;
	vx_t = 0.;
	vy_t = 0.;
	Sxx_t = 0.;
	Sxy_t = 0.;
	Syy_t = 0.;
}

void particle::copy_into(particle *p) const {
	p->idx = idx;
	p->hash = hash;

	p->m = m;
	p->rho = rho;
	p->h = h;
	p->quad_weight = quad_weight;
	p->p = this->p;

	p->x = x;
	p->y = y;

	p->X = X;
	p->Y = Y;

	p->vx = vx;
	p->vy = vy;

	p->Sxx = Sxx;
	p->Sxy = Sxy;
	p->Syy = Syy;

	p->fcx = fcx;
	p->fcy = fcy;

	p->num_nbh = num_nbh;
	memcpy(p->nbh, nbh, sizeof(unsigned int)*num_nbh);
	memcpy(p->w,   w,   sizeof(kernel_result)*num_nbh);
}

//-----------------------------------------------------

particle_ul::particle_ul() : particle() {
}

particle_ul::particle_ul(unsigned int idx) : particle(idx){
}

particle_ul::~particle_ul() {};

void particle_ul::reset() {
	particle::reset();
}
void particle_ul::copy_into(particle *p) const {
	particle::copy_into(p);
}

//-----------------------------------------------------

particle_tl::particle_tl() : particle() {
}

particle_tl::particle_tl(unsigned int idx) : particle(idx){
}

particle_tl::~particle_tl() {};

void particle_tl::reset() {
	particle::reset();
}
void particle_tl::copy_into(particle *p) const {
	particle::copy_into(p);

	assert(static_cast<particle_tl*>(p) != NULL);
	particle_tl *p_other = static_cast<particle_tl*>(p);

	p_other->Pxx = this->Pxx;
	p_other->Pxy = this->Pxy;
	p_other->Pyx = this->Pyx;
	p_other->Pyy = this->Pyy;
}

//-----------------------------------------------------

particle_corotational::particle_corotational() : particle_tl() {
}

particle_corotational::particle_corotational(unsigned int idx) : particle_tl(idx){
}

particle_corotational::~particle_corotational() {};

void particle_corotational::reset() {
	particle_tl::reset();
}
void particle_corotational::copy_into(particle *p) const {
	particle_tl::copy_into(p);

	assert(static_cast<particle_corotational*>(p) != NULL);
	particle_corotational *p_other = static_cast<particle_corotational*>(p);

	p_other->Rxx = this->Rxx;
	p_other->Rxy = this->Rxy;
	p_other->Ryx = this->Ryx;
	p_other->Ryy = this->Ryy;
}

//-----------------------------------------------------

particle_potential::particle_potential() : particle_tl() {
}

particle_potential::particle_potential(unsigned int idx) : particle_tl(idx){
}

particle_potential::~particle_potential() {};

void particle_potential::reset() {
	particle_tl::reset();
}
void particle_potential::copy_into(particle *p) const {
	particle_tl::copy_into(p);

	assert(static_cast<particle_potential*>(p) != NULL);
	particle_potential *p_other = static_cast<particle_potential*>(p);

	memcpy(p_other->wji, wji, sizeof(kernel_result)*num_nbh);
}

//-----------------------------------------------------

particle_tl_weak::particle_tl_weak() : particle_tl() {
}

particle_tl_weak::particle_tl_weak(unsigned int idx) : particle_tl(idx){
}

particle_tl_weak::~particle_tl_weak() {};

void particle_tl_weak::reset() {
	particle_tl::reset();
}
void particle_tl_weak::copy_into(particle *p) const {
	particle_tl::copy_into(p);
}
//-----------------------------------------------------

particle_ul_weak::particle_ul_weak() : particle_ul() {
}

particle_ul_weak::particle_ul_weak(unsigned int idx) : particle_ul(idx) {
}

particle_ul_weak::~particle_ul_weak() {};

void particle_ul_weak::reset() {
	particle_ul::reset();
}

void particle_ul_weak::copy_into(particle *p) const  {
	particle_ul::copy_into(p);

	assert(static_cast<particle_ul_weak*>(p) != NULL);
	particle_ul_weak *p_other = static_cast<particle_ul_weak*>(p);

	p_other->num_cand_nbh = num_cand_nbh;
	memcpy(p_other->cand_nbh, cand_nbh, num_cand_nbh*sizeof(unsigned int));
}

//-----------------------------------------------------

particle_monaghan::particle_monaghan() : particle_ul() {
}

particle_monaghan::particle_monaghan(unsigned int idx) : particle_ul(idx) {
}

particle_monaghan::~particle_monaghan() {};

void particle_monaghan::reset() {
	particle_ul::reset();
}

void particle_monaghan::copy_into(particle *p) const  {
	particle_ul::copy_into(p);

	assert(static_cast<particle_monaghan*>(p) != NULL);
	particle_monaghan *p_other = static_cast<particle_monaghan*>(p);

	p_other->Rxx = this->Rxx;
	p_other->Rxy = this->Rxy;
	p_other->Ryy = this->Ryy;
}
