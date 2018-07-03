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

#ifndef ROTATION_H_
#define ROTATION_H_

#include "../particle.h"
#include "../body.h"

#include "../correctors.h"
#include "../cont_mech.h"
#include "../material.h"
#include "../derivatives.h"
#include "../godunov.h"
#include "../riemann_integration.h"
#include "../convex_hull.h"
#include "../boundary.h"
#include "../interpolation.h"
#include "../corotational.h"
#include "../potential.h"

#include "neighbor_search.h"
#include "utils.h"

std::vector<body<particle_monaghan>>     rotation_monaghan(unsigned int nbox);		// wobbles (expected)
std::vector<body<particle_ul>>           rotation_godunov(unsigned int nbox);		// wobbles (expected)
std::vector<body<particle_corotational>> rotation_corot(unsigned int nbox);			// works
std::vector<body<particle_potential>>    rotation_potential(unsigned int nbox);		// works
std::vector<body<particle_tl>>           rotation_tl(unsigned int nbox);			// fails (expected? LSDYNA fails too)
std::vector<body<particle_ul_weak>>      rotation_ul_fem(unsigned int nbox);		// works
std::vector<body<particle_tl_weak>>      rotation_tl_fem(unsigned int nbox);		// works
std::vector<body<particle_ul_weak>> 	 rotation_ul_weak(unsigned int nbox);		// works
std::vector<body<particle_tl_weak>>      rotation_tl_weak(unsigned int nbox);		// works

#endif /* ROTATION_H_ */
