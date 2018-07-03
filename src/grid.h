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

//standard cell list implementation. is used to turn cell lists into verlet lists which are then used in turn
// in derivatives.h

#ifndef GRID_H_
#define GRID_H_

#include "particle.h"

#include <vector>
#include <float.h>
#include <math.h>
#include <assert.h>

template <typename T>
class grid {

public:
	void assign_hashes(T ** particles, unsigned int n) const {
		for (unsigned int i = 0; i < n; i++) {
			unsigned int ix = (unsigned int) ((particles[i]->x - m_bbmin_x)/m_dx);
			unsigned int iy = (unsigned int) ((particles[i]->y - m_bbmin_y)/m_dx);
			particles[i]->hash = ix*m_ny + iy;
		}
	}

	void unhash(int idx, unsigned int &i, unsigned int &j) const {
		i = idx/(m_ny);
		j = idx-(i)*(m_ny);
	}

	std::vector<int> get_cells( const T * const * const particles, unsigned int n) const {

		//needs to be sorted for this to work
		for (unsigned int i = 0; i < n-1; i++) {
			assert(particles[i]->hash <= particles[i+1]->hash);
		}

		std::vector<int> cells(m_num_cell+1);
		std::fill(cells.begin(), cells.end(), -1);

		cells[particles[0]->hash] = 0;
		for (unsigned int i = 0; i < n-1; i++) {
			if (particles[i]->hash != particles[i+1]->hash) {
				cells[particles[i+1]->hash] = i+1;
			}
		}
		cells[particles[n-1]->hash+1] = n;

		//empty boxes are now set to -1
		//in order to iterate through a cell by [cells(cell_index),...,cells(cell_index+1)[
		//those need to be fixed by propagating a "fix" value from the right
		//(such that the above range will just be empty)

		unsigned int fix = n;
		for (auto it = cells.rbegin(); it != cells.rend(); ++it) {
			if (*it == -1) {
				*it = fix;
			} else {
				fix = *it;
			}
		}

		return cells;
	}

	void update_geometry(const T * const * const particles, unsigned int n, double kernel_width) {
		m_bbmin_x = DBL_MAX; m_bbmax_x = -DBL_MAX;
		m_bbmin_y = DBL_MAX; m_bbmax_y = -DBL_MAX;

		double h_max = -DBL_MAX;

		for (unsigned int i = 0; i < n; i++) {
			m_bbmin_x = fmin(particles[i]->x, m_bbmin_x);
			m_bbmin_y = fmin(particles[i]->y, m_bbmin_y);
			m_bbmax_x = fmax(particles[i]->x, m_bbmax_x);
			m_bbmax_y = fmax(particles[i]->y, m_bbmax_y);

			h_max = fmax(particles[i]->h, h_max);
		}

		//some nudging to prevent round off errors

		m_bbmin_x -= 1e-6;
		m_bbmin_y -= 1e-6;
		m_bbmax_x += 1e-6;
		m_bbmax_y += 1e-6;

		m_dx = h_max*kernel_width;

		m_lx = m_bbmax_x - m_bbmin_x;
		m_ly = m_bbmax_y - m_bbmin_y;

		m_nx = ceil(m_lx/m_dx);
		m_ny = ceil(m_ly/m_dx);
		m_num_cell = m_nx*m_ny;

		assert(m_num_cell < n);
	}

	void construct_verlet_lists(T **particles, unsigned int n, double kernel_width) {
		std::vector<int> cells = get_cells(particles, n);

		for (unsigned int b = 0; b < m_num_cell; b++) {
			unsigned int gi = 0; unsigned int gj = 0;

			unhash((int) b, gi, gj);

			int low_i  = (int) gi-1 < 0 ? 0 : (int) gi-1;
			int low_j  = (int) gj-1 < 0 ? 0 : (int) gj-1;
			int high_i = gi+2 > m_nx ? m_nx : gi+2;
			int high_j = gj+2 > m_ny ? m_ny : gj+2;

			for (int i = cells[b]; i < cells[b+1]; i++) {

				particle *p1 = particles[i];
				unsigned int nbh_iter = 0;

				double hi = p1->h;

				double xi = p1->x;
				double yi = p1->y;

				double radius2 = hi*hi*2*2;

				for (int ni = low_i; ni < high_i; ni++) {
					for (int nj = low_j; nj < high_j; nj++) {

						for (int j = cells[ni*m_ny+nj]; j < cells[ni*m_ny+nj+1]; j++) {

							particle *p2 = particles[j];

							double xj = p2->x;
							double yj = p2->y;

							double xij = xi-xj;
							double yij = yi-yj;

							double r2 = xij*xij + yij*yij;

							if (r2 <= radius2 + 1e-6) {
								particles[i]->nbh[nbh_iter] = j;
								nbh_iter++;
							}

						}

					}
				}

				assert(nbh_iter < MAX_NBH);
				assert(nbh_iter > 0);

				if (nbh_iter == 0) {
					printf("alarm, particle with no neighbors found!\n");
				}

				particles[i]->num_nbh = nbh_iter;
			}
		}
	}

private:
	double m_dx;
	double m_lx, m_ly;
	unsigned int m_nx, m_ny;
	unsigned int m_num_cell;

	double m_bbmin_x, m_bbmax_x;
	double m_bbmin_y, m_bbmax_y;

	std::vector<int> m_cells;
};

#endif /* GRID_H_ */
