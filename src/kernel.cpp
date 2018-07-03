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

#include "kernel.h"

kernel_result cubic_spline(double xi, double yi, double xj, double yj, double h) {
	double h1 = 1./h;
	double xij = xi-xj;
	double yij = yi-yj;
	double rij = sqrt(xij*xij + yij*yij);

	double fac = 10*(M_1_PI)/7.0*h1*h1;
	double q = rij*h1;

	kernel_result w;
	w.w = 0.;
	w.w_x = 0.;
	w.w_y = 0.;

	if (q >= 2.) {
		return w;
	} else if (q >= 1.) {
		w.w =  fac*(0.25*(2-q)*(2-q)*(2-q));
		if (rij > 1e-12) {
			double der = -0.75*(2-q)*(2-q) * h1/rij;
			w.w_x = xij*der*fac;
			w.w_y = yij*der*fac;
		}
		return w;
	} else {
		w.w = fac*(1 - 1.5*q*q*(1-0.5*q));
		if (rij > 1e-12) {
			double der = -3.0*q*(1-0.75*q) * h1/rij;
			w.w_x = xij*der*fac;
			w.w_y = yij*der*fac;
		}
		return w;
	}

	return w;
}
