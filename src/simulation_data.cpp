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

#include "simulation_data.h"

physical_constants::physical_constants(double nu, double E, double rho0) : m_nu(nu), m_E(E), m_rho0(rho0) {}

physical_constants::physical_constants() {}

double physical_constants::nu() const {
	return m_nu;
}

double physical_constants::E() const {
	return m_E;
}

double physical_constants::G() const {
	return m_E/(2.*(1.+m_nu));
}

double physical_constants::K()  const {
	double G = this->G();
	return 2.0*G*(1+m_nu)/(3*(1-2*m_nu));
}

double physical_constants::rho0() const {
	return m_rho0;
}

double physical_constants::c0() const {
	return sqrt(K()/m_rho0);
}

//---------------------------------------------------------------------------

double constants_monaghan::mghn_wdeltap() const {
	return m_mghn_wdeltap;
}

double constants_monaghan::mghn_corr_exp() const {
	return m_mghn_corr_exp;
}

double constants_monaghan::mghn_eps() const {
	return m_mghn_eps;
}

constants_monaghan::constants_monaghan(double wdeltap, double corr_exp, double eps) :
		m_mghn_wdeltap(wdeltap),
		m_mghn_corr_exp(corr_exp),
		m_mghn_eps(eps) {}

constants_monaghan::constants_monaghan() {}

//---------------------------------------------------------------------------

double constants_artificial_viscosity::artvisc_alpha() const {
	return m_artvisc_alpha;
}

double constants_artificial_viscosity::artvisc_beta() const {
	return m_artvisc_beta;
}

double constants_artificial_viscosity::artvisc_eta() const {
	return m_artvisc_eta;
}

constants_artificial_viscosity::constants_artificial_viscosity(double alpha, double beta, double eta) :
	m_artvisc_alpha(alpha),
	m_artvisc_beta(beta),
	m_artvisc_eta(eta) {}

constants_artificial_viscosity::constants_artificial_viscosity() {}

//---------------------------------------------------------------------------

correction_constants::correction_constants(constants_monaghan monaghan_constants, constants_artificial_viscosity constants_art_visc, double xsph_eps) :
		m_xsph_eps(xsph_eps),
		m_constants_monaghan(monaghan_constants),
		m_constants_art_visc(constants_art_visc)
		 {}

correction_constants::correction_constants(constants_monaghan monaghan_constants, constants_artificial_viscosity constants_art_visc, double xsph_eps, bool smooth) :
		m_xsph_eps(xsph_eps),
		m_do_smooth(smooth),
		m_constants_monaghan(monaghan_constants),
		m_constants_art_visc(constants_art_visc)
		{}

correction_constants::correction_constants() {}


double correction_constants::xsph_eps() const {
	return m_xsph_eps;
}

bool correction_constants::do_smooth() const {
	return m_do_smooth;
}

constants_monaghan correction_constants::get_monaghan_const() const {
	return m_constants_monaghan;
}

constants_artificial_viscosity correction_constants::get_art_visc_const() const {
	return m_constants_art_visc;
}

//---------------------------------------------------------------------------

simulation_data::simulation_data(physical_constants physical_constants, correction_constants correction_constants)
	: m_physical_constants(physical_constants), m_correction_constants(correction_constants) {};

simulation_data::simulation_data() {}

physical_constants simulation_data::get_physical_constants() const {
	return m_physical_constants;
}

correction_constants simulation_data::get_correction_constants() const {
	return m_correction_constants;
}
