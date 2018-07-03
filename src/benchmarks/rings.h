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

//rubber ring impact according to
//Gray, J. P., J. J. Monaghan, and R. P. Swift. "SPH elastic dynamics." Computer methods in applied mechanics and engineering 190.49-50 (2001): 6641-6662.

#ifndef RINGS_H_
#define RINGS_H_

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

// not using contact algorithm
std::vector<body<particle_monaghan>> rings_monaghan(unsigned int nbox);	//works!
std::vector<body<particle_ul>>       rings_godunov(unsigned int nbox);  //works!
// using contact algorithm
std::vector<body<particle_monaghan>>     rings_monaghan_contact(unsigned int nbox);  // works (needed velocity decrease)
std::vector<body<particle_ul>>           rings_godunov_contact(unsigned int nbox);   // works
std::vector<body<particle_corotational>> rings_corot_contact(unsigned int nbox);     // works
std::vector<body<particle_potential>>    rings_potential_contact(unsigned int nbox); // ridicilously small timestep needed? can not be improved by rate model
std::vector<body<particle_tl>>           rings_tl_contact(unsigned int nbox);        // works but needs high amount of artificial viscosity
std::vector<body<particle_ul_weak>>      rings_ul_fem_contact(unsigned int nbox);	 // works
std::vector<body<particle_tl_weak>>      rings_tl_fem_contact(unsigned int nbox);    // works, but only with rate material model
std::vector<body<particle_ul_weak>> 	 rings_ul_weak_contact(unsigned int nbox);	 // works, kinda dissipative
std::vector<body<particle_tl_weak>>      rings_tl_weak_contact(unsigned int nbox);	 // works  but only with rate material model (?)

#endif /* RINGS_H_ */
