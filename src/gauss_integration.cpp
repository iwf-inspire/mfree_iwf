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

#include "gauss_integration.h"

void gauss_int::four_node_quadrilateral(double xi, double eta, glm::dvec4 &N, glm::dmat2x4 &N_xi) {
	double X1 = -1; double Y1 = -1;
	double X2 =  1; double Y2 = -1;
	double X3 =  1; double Y3 =  1;
	double X4 = -1; double Y4 =  1;


	N[0] = 0.25*(1 + X1*xi)*(1 + Y1*eta);
	N[1] = 0.25*(1 + X2*xi)*(1 + Y2*eta);
	N[2] = 0.25*(1 + X3*xi)*(1 + Y3*eta);
	N[3] = 0.25*(1 + X4*xi)*(1 + Y4*eta);

	N_xi[0][0] = 0.25*X1*(1 + Y1*eta); N_xi[1][0] = 0.25*Y1*(1 + X1*xi);
	N_xi[0][1] = 0.25*X2*(1 + Y2*eta); N_xi[1][1] = 0.25*Y2*(1 + X2*xi);
	N_xi[0][2] = 0.25*X3*(1 + Y3*eta); N_xi[1][2] = 0.25*Y3*(1 + X3*xi);
	N_xi[0][3] = 0.25*X4*(1 + Y4*eta); N_xi[1][3] = 0.25*Y4*(1 + X4*xi);
}

void gauss_int::get_weight_pos_1d(std::vector<double> &p1, std::vector<double> &w1, unsigned int order) {
	assert(order > 0);
	assert(order < 9);

	p1.resize(order);
	w1.resize(order);

	switch (order) {
	case 1:
		p1[0] = 0.;
		w1[0] = 2.;
		break;
	case 2:
		p1[0] = 0.577350269189626;
		p1[1] =-0.577350269189626;

		w1[0] = 1.000000000000000;
		w1[1] = 1.000000000000000;
		break;
	case 3:
		p1[0] = 0.774596669241483;
		p1[1] =-0.774596669241483;
		p1[2] = 0.000000000000000;

		w1[0] = 0.555555555555556;
		w1[1] = 0.555555555555556;
		w1[2] = 0.888888888888889;
		break;
	case 4:
		p1[0] = 0.861134311594053;
		p1[1] =-0.861134311594053;
		p1[2] = 0.339981043584856;
		p1[3] =-0.339981043584856;

		w1[0] = 0.347854845137454;
		w1[1] = 0.347854845137454;
		w1[2] = 0.652145154862546;
		w1[3] = 0.652145154862546;
		break;
	case 5:
		p1[0] = 0.906179845938664;
		p1[1] =-0.906179845938664;
		p1[2] = 0.538469310105683;
		p1[3] =-0.538469310105683;
		p1[4] = 0.000000000000000;

		w1[0] = 0.236926885056189;
		w1[1] = 0.236926885056189;
		w1[2] = 0.478628670499366;
		w1[3] = 0.478628670499366;
		w1[4] = 0.568888888888889;
		break;
	case 6:
		p1[0] = 0.932469514203152;
		p1[1] =-0.932469514203152;
		p1[2] = 0.661209386466265;
		p1[3] =-0.661209386466265;
		p1[4] = 0.238619186003152;
		p1[5] =-0.238619186003152;

		w1[0] = 0.171324492379170;
		w1[1] = 0.171324492379170;
		w1[2] = 0.360761573048139;
		w1[3] = 0.360761573048139;
		w1[4] = 0.467913934572691;
		w1[5] = 0.467913934572691;
		break;
	case 7:
		p1[0] =  0.949107912342759;
		p1[1] = -0.949107912342759;
		p1[2] =  0.741531185599394;
		p1[3] = -0.741531185599394;
		p1[4] =  0.405845151377397;
		p1[5] = -0.405845151377397;
		p1[6] =  0.000000000000000;

		w1[0] = 0.129484966168870;
		w1[1] = 0.129484966168870;
		w1[2] = 0.279705391489277;
		w1[3] = 0.279705391489277;
		w1[4] = 0.381830050505119;
		w1[5] = 0.381830050505119;
		w1[6] = 0.417959183673469;
		break;
	case 8:
		p1[0] =  0.960289856497536;
		p1[1] = -0.960289856497536;
		p1[2] =  0.796666477413627;
		p1[3] = -0.796666477413627;
		p1[4] =  0.525532409916329;
		p1[5] = -0.525532409916329;
		p1[6] =  0.183434642495650;
		p1[7] = -0.183434642495650;

		w1[0] = 0.101228536290376;
		w1[1] = 0.101228536290376;
		w1[2] = 0.222381034453374;
		w1[3] = 0.222381034453374;
		w1[4] = 0.313706645877887;
		w1[5] = 0.313706645877887;
		w1[6] = 0.362683783378362;
		w1[7] = 0.362683783378362;
		break;
	}
}

void gauss_int::get_weight_pos_2d(std::vector<std::tuple<double, double> > &pt, std::vector<double> &w, unsigned int order) {

	std::vector<double> p1;
	std::vector<double> w1;

	get_weight_pos_1d(p1, w1, order);

	for (unsigned int i = 0; i < order; i++) {
		for (unsigned int j = 0; j < order; j++) {
			pt.push_back(std::tuple<double, double>(p1[i],p1[j]));
			w.push_back(w1[i]*w1[j]);
		}
	}
}

gauss_int::gauss_int(std::vector<element> elements, unsigned int order) : m_elements(elements), m_order(order) {
		m_local_gp_pos_x  = (double*) calloc(order*order, sizeof(double));
		m_local_gp_pos_y  = (double*) calloc(order*order, sizeof(double));
		m_local_gp_weight = (double*) calloc(order*order, sizeof(double));

		std::vector<std::tuple<double, double> > pt;
		std::vector<double> w;
		get_weight_pos_2d(pt, w, order);

		for (unsigned int i = 0; i < order*order; i++) {
			m_local_gp_pos_x[i]  = std::get<0>(pt[i]);
			m_local_gp_pos_y[i]  = std::get<1>(pt[i]);
			m_local_gp_weight[i] = w[i];
		}

	}
