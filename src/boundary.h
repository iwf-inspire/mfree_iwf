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

//boundary data structure and methods for contact algorithm
//stores normals as circular tour. normals are updated according to simple geometric principles, see figure (6)
//and formulas 164-168 in the paper

#ifndef BOUNDARY_H_
#define BOUNDARY_H_

#include "kd_tree.h"
#include <glm/glm.hpp>

template<typename T>
class boundary {
private:
	unsigned int *m_indexes = 0;
	kd_point *m_boundary = 0;
	unsigned int m_num_boundary = 0;

public:
	boundary() {}

	boundary(std::vector<unsigned int> indexes, T** particles) {
		m_num_boundary = indexes.size();
		m_boundary = new kd_point[m_num_boundary];
		m_indexes = new unsigned int[m_num_boundary];

		std::copy(indexes.begin(), indexes.end(), m_indexes);
		for (unsigned int i = 0; i < m_num_boundary; i++) {
			m_boundary[i].x = particles[i]->x;
			m_boundary[i].y = particles[i]->y;
		}
	}

	void update_normals() {
		for (unsigned int i = 0; i < m_num_boundary; i++) {
			int left = ((int) i) - 1;
			int right = i+1;

			if (left < 0)
				left = m_num_boundary-1;

			if (right >= (int) m_num_boundary)
				right = 0;

			glm::dvec2 left_p(m_boundary[left].x, m_boundary[left].y);
			glm::dvec2 p(m_boundary[i].x, m_boundary[i].y);
			glm::dvec2 right_p(m_boundary[right].x, m_boundary[right].y);

			glm::dvec2 tl = left_p - p;
			glm::dvec2 tr = right_p - p;

			glm::dvec2 nl = glm::normalize(glm::dvec2(-tl.y,  tl.x));
			glm::dvec2 nr = glm::normalize(glm::dvec2( tr.y, -tr.x));

			glm::dvec2 n = 0.5*(nl + nr);

			m_boundary[i].nx = n.x;
			m_boundary[i].ny = n.y;
		}
	}

	kd_point *points() const {
		return m_boundary;
	}

	const unsigned int *indexes() const {
		return (const unsigned int *) m_indexes;
	}

	unsigned int num_bound() const {
		return m_num_boundary;
	}
};

#endif /* BOUNDARY_H_ */
