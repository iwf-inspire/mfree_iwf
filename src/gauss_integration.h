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

//perform gauss quadrature for FEM and the (meshfree) weak TL method
//	adapted from the companian code to Nguyen, Vinh Phu, et al. "Meshless methods: a review and computer implementation aspects." Mathematics and computers in simulation 79.3 (2008): 763-813.

#ifndef GAUSS_INTEGRATION_H_
#define GAUSS_INTEGRATION_H_

#include "particle.h"
#include "element.h"

#include <vector>
#include <tuple>
#include <assert.h>
#include <glm/glm.hpp>

class gauss_int {
private:
	std::vector<element> m_elements;
	double *m_local_gp_pos_x;
	double *m_local_gp_pos_y;
	double *m_local_gp_weight;
	unsigned int m_order;

	static void four_node_quadrilateral(double xi, double eta, glm::dvec4 &N, glm::dmat2x4 &N_xi);
	static void get_weight_pos_1d(std::vector<double> &p1, std::vector<double> &w1, unsigned int order);
	static void get_weight_pos_2d(std::vector<std::tuple<double, double> > &pt, std::vector<double> &w, unsigned int order);

public:
	gauss_int(std::vector<element> elements, unsigned int order);

	template<typename T>
	void update_gauss_points(T** particles, T** gauss_points) {
		unsigned int global_gp_iter = 0;
		unsigned int num_gp_2d = m_order*m_order;

		for (auto it = m_elements.begin(); it != m_elements.end(); ++it) {
			unsigned int node_idx1 = it->ll_idx;
			unsigned int node_idx2 = it->lr_idx;
			unsigned int node_idx3 = it->tr_idx;
			unsigned int node_idx4 = it->tl_idx;

			unsigned int node_idx[4] = {node_idx1, node_idx2, node_idx3, node_idx4};

			glm::dmat2x4 m_node;
			m_node[0][0] = particles[node_idx1]->x; m_node[1][0] = particles[node_idx1]->y;
			m_node[0][1] = particles[node_idx2]->x; m_node[1][1] = particles[node_idx2]->y;
			m_node[0][2] = particles[node_idx3]->x; m_node[1][2] = particles[node_idx3]->y;
			m_node[0][3] = particles[node_idx4]->x; m_node[1][3] = particles[node_idx4]->y;

			glm::dmat2x4 m_node_par;
			m_node_par[0][0] = -1; m_node_par[1][0] = -1;
			m_node_par[0][1] = +1; m_node_par[1][1] = -1;
			m_node_par[0][2] = +1; m_node_par[1][2] = +1;
			m_node_par[0][3] = -1; m_node_par[1][3] = +1;

			for (unsigned int i = 0; i < num_gp_2d; i++) {
				double local_x = m_local_gp_pos_x[i];
				double local_y = m_local_gp_pos_y[i];
				double local_w = m_local_gp_weight[i];

				glm::dvec4 N;
				glm::dmat2x4 N_xi;
				four_node_quadrilateral(local_x, local_y, N, N_xi);

				//position
				glm::dvec2 global_pos = glm::transpose(m_node)*N;

				gauss_points[global_gp_iter]->x = global_pos[0];
				gauss_points[global_gp_iter]->y = global_pos[1];

				gauss_points[global_gp_iter]->X = global_pos[0];
				gauss_points[global_gp_iter]->Y = global_pos[1];

				glm::dmat2x2 x_xi(0.);
				for (unsigned int I = 0; I < 4; I++) {
		            double x_xi11 = 0.25*m_node[0][I]*m_node_par[0][I]*(1+m_node_par[1][I]*local_y);
		            double x_xi12 = 0.25*m_node[0][I]*m_node_par[1][I]*(1+m_node_par[0][I]*local_x);
		            double x_xi21 = 0.25*m_node[1][I]*m_node_par[0][I]*(1+m_node_par[1][I]*local_y);
		            double x_xi22 = 0.25*m_node[1][I]*m_node_par[1][I]*(1+m_node_par[0][I]*local_x);

		            x_xi += glm::dmat2x2(x_xi11, x_xi21, x_xi12, x_xi22);
				}

				//weight
				double J = -glm::determinant(x_xi);
				double global_gp_weigth = local_w*J;
				assert(global_gp_weigth > 0.);

				gauss_points[global_gp_iter]->quad_weight  = global_gp_weigth;

				//function values
				glm::dmat2x2 x_xi_inv = glm::inverse(x_xi);
				glm::dmat2x4 N_x = N_xi*x_xi_inv;

				gauss_points[global_gp_iter]->num_nbh = 4;
				for (unsigned int I = 0; I < 4; I++) {
					gauss_points[global_gp_iter]->nbh[I] = node_idx[I];
					gauss_points[global_gp_iter]->w[I].w = N[I];

					gauss_points[global_gp_iter]->w[I].w_x = N_x[0][I];
					gauss_points[global_gp_iter]->w[I].w_y = N_x[1][I];
				}

				global_gp_iter++;
			}
		}
	}
};


#endif /* GAUSS_INTEGRATION_H_ */
