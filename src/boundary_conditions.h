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

// methods to apply (Dirichlet) boundary conditions. takes a tuple of indices into the particle array and a function pointer
// that may, in principle, apply any transformation of the state and time derived variables of the particles at the indices

#ifndef BOUNDARY_CONDITIONS_H_
#define BOUNDARY_CONDITIONS_H_

#include <vector>
#include <tuple>

template<typename T>
class boundary_conditions {
public:
	boundary_conditions() {};

	void add_boundary_condition(const std::vector<unsigned int> boundary_indexes, void(*boundary_fun)(T*)) {
		m_boundary_conditions.push_back(std::tuple<std::vector<unsigned int>, void(*)(T*)>(boundary_indexes, boundary_fun));
	}

	void carry_out_boundary_conditions(T **particles, unsigned int np) const {
		for (auto it = m_boundary_conditions.begin(); it != m_boundary_conditions.end(); ++it) {
			std::vector<unsigned int> cur_indices = std::get<0>(*it);
			auto cur_bc = std::get<1>(*it);

			for (auto jt = cur_indices.begin(); jt != cur_indices.end(); ++jt) {
				unsigned int idx = *jt;
				cur_bc(particles[idx]);
			}
		}
	}

private:
	std::vector<std::tuple <std::vector<unsigned int>, void(*)(T*)> > m_boundary_conditions;

};

#endif /* BOUNDARY_CONDITIONS_H_ */
