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

//basic control flow is like this
//	benchmark files construct some body or bodies
//	body and time stepper are friends
//	body "knows" how to compute time derivatives
//	time stepper has body compute time derivatives as needed by the time stepping scheme (4 times in case of runge kutta 4)


#include <iostream>
#include <stdlib.h>
#include <fenv.h>

#include "body.h"
#include "particle.h"
#include "contact.h"

#include "vtk_writer.h"

#include "benchmarks/rings.h"
#include "benchmarks/tensile.h"
#include "benchmarks/rotation.h"

int main() {
	feenableexcept(FE_INVALID | FE_OVERFLOW);
	system("rm ./results/*.vtk");
	system("rm ./results/*.fbg*");

//	std::vector<body<particle_ul_weak>>      bodies = rings_ul_fem_contact(80);
	std::vector<body<particle_monaghan>>     bodies = rings_monaghan(80);
//	std::vector<body<particle_monaghan>>     bodies = rings_monaghan_contact(80);
//	std::vector<body<particle_ul>>           bodies = rings_godunov(80);
//	std::vector<body<particle_ul>>           bodies = rings_godunov_contact(80);
//	std::vector<body<particle_corotational>> bodies = rings_corot_contact(80);
//	std::vector<body<particle_potential>>    bodies = rings_potential_contact(80);
//	std::vector<body<particle_tl>>           bodies = rings_tl_contact(80);
//	std::vector<body<particle_tl_weak>>      bodies = rings_tl_fem_contact(80);
//	std::vector<body<particle_ul_weak>>      bodies = rings_ul_weak_contact(80);
//	std::vector<body<particle_tl_weak>>      bodies = rings_tl_weak_contact(80);

//	std::vector<body<particle_ul_weak>>      bodies = rotation_ul_fem(20);
//	std::vector<body<particle_monaghan>>     bodies = rotation_monaghan(20);
//	std::vector<body<particle_ul>>           bodies = rotation_godunov(20);
//	std::vector<body<particle_corotational>> bodies = rotation_corot(20);
//	std::vector<body<particle_potential>>    bodies = rotation_potential(20);
//	std::vector<body<particle_tl>>           bodies = rotation_tl(20);
//	std::vector<body<particle_tl_weak>>      bodies = rotation_tl_fem(20);
//	std::vector<body<particle_ul_weak>>      bodies = rotation_ul_weak(20);
//	std::vector<body<particle_tl_weak>>      bodies = rotation_tl_weak(20);

//	std::vector<body<particle_ul_weak>>      bodies = tensile_ul_fem(21);
//	std::vector<body<particle_monaghan>>     bodies = tensile_monaghan(21);
//	std::vector<body<particle_monaghan>>     bodies = tensile_monaghan_recreate();
//	std::vector<body<particle_ul>>           bodies = tensile_godunov(21);
//	std::vector<body<particle_corotational>> bodies = tensile_corot(21);
//	std::vector<body<particle_potential>>    bodies = tensile_potential(21);
//	std::vector<body<particle_tl>>           bodies = tensile_tl(21);
//	std::vector<body<particle_tl_weak>>      bodies = tensile_tl_fem(21);
//	std::vector<body<particle_ul_weak>>      bodies = tensile_ul_weak(21);
//	std::vector<body<particle_tl_weak>>      bodies = tensile_tl_weak(21);

	unsigned int num_bodies = bodies.size();

	simulation_time *time = &simulation_time::getInstance();

	unsigned int step =      0;
	unsigned int freq =     10;

	while(true) {

		for (unsigned int i = 0; i < num_bodies; i++) {
			bodies[i].step();
			bodies[i].update_boundary();
		}

		for (unsigned int i = 0; i < num_bodies; i++) {
			for (unsigned int j = 0; j < num_bodies; j++) {
				if (i != j) {
					contact_apply_body_to_boundary(bodies[i], bodies[j]);
				}
			}
		}

		if (step % freq == 0) {
			vtk_writer_write(bodies, step);
		}

		printf("%d: %f\n", step, time->get_time());
		step++;
	}

	return EXIT_SUCCESS;
}
