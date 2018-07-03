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

//structs that hold all immutable simulation data (like physcial constants)

#ifndef SIMULATION_DATA_H_
#define SIMULATION_DATA_H_

#include <math.h>

class physical_constants {
private:
	double m_nu = 0.;
	double m_E = 0.;
	double m_rho0 = 0.;

public:
	physical_constants(double nu, double E, double rho0);
	physical_constants();
	double nu() const;
	double E() const;
	double G() const;
	double K() const;
	double rho0() const;
	double c0() const;
};

class constants_monaghan {
	double m_mghn_wdeltap = 0.;
	double m_mghn_corr_exp = 0.;
	double m_mghn_eps = 0.;

public:
	double mghn_wdeltap() const;
	double mghn_corr_exp() const;
	double mghn_eps() const;

	constants_monaghan(double wdeltap, double corr_exp, double eps);
	constants_monaghan();
};

class constants_artificial_viscosity {
	double m_artvisc_alpha = 0.;
	double m_artvisc_beta = 0.;
	double m_artvisc_eta = 0.;

public:
	double artvisc_alpha() const;
	double artvisc_beta() const;
	double artvisc_eta() const;

	constants_artificial_viscosity(double alpha, double beta, double eta);
	constants_artificial_viscosity();
};

class correction_constants {

private:
	double m_xsph_eps = 0.;
	bool m_do_smooth = false;
	constants_monaghan m_constants_monaghan;
	constants_artificial_viscosity m_constants_art_visc;

public:
	correction_constants(constants_monaghan monaghan_constants, constants_artificial_viscosity constants_art_visc, double xsph_eps);
	correction_constants(constants_monaghan monaghan_constants, constants_artificial_viscosity constants_art_visc, double xsph_eps, bool smooth);
	correction_constants();

	double xsph_eps() const;
	bool do_smooth() const;
	constants_monaghan get_monaghan_const() const;
	constants_artificial_viscosity get_art_visc_const() const;
};

class simulation_data {
private:
	physical_constants   m_physical_constants;
	correction_constants m_correction_constants;

public:
	simulation_data();
	simulation_data(physical_constants physical_constants, correction_constants correction_constants);
	physical_constants   get_physical_constants() const;
	correction_constants get_correction_constants() const;
};

#endif /* SIMULATION_DATA_H_ */
