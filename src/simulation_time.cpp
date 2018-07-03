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

#include "simulation_time.h"

simulation_time& simulation_time::getInstance() {
    static simulation_time instance;
    return instance;
}

double simulation_time::get_time() const {
	return m_time;
}

double simulation_time::get_dt() const {
	return m_dt;
}

 simulation_time::simulation_time() : m_time(0.), m_dt(0.) {
 }

 void simulation_time::increment_time() {
	 m_time += m_dt;
 }


