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

//quick and dirty writer for xml based legacy vtk format. use paraview to see simulation results

#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include "particle.h"

#include <stdio.h>

template <typename T>
void vtk_writer_write(const std::vector<body<T>> &bodies, unsigned int step) {

	unsigned int total_num_part = 0;
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		total_num_part += it->get_num_part();
	}

	char buf[256];
	sprintf(buf, "results/out_%06d.vtk", step);
	FILE *fp = fopen(buf, "w+");

	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "mfree iwf\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "\n");

	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", total_num_part);
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f %f %f\n", particles[i]->x, particles[i]->y, 0.);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELLS %d %d\n", total_num_part, 2*total_num_part);
	for (unsigned int i = 0; i < total_num_part; i++) {
		fprintf(fp, "%d %d\n", 1, i);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", total_num_part);
	for (unsigned int i = 0; i < total_num_part; i++) {
		fprintf(fp, "%d\n", 1);
	}
	fprintf(fp, "\n");

	fprintf(fp, "POINT_DATA %d\n", total_num_part);
	fprintf(fp, "SCALARS h float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f\n", particles[i]->h);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS density float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f\n", particles[i]->rho);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS Sxx float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f\n", particles[i]->Sxx);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS Sxy float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f\n", particles[i]->Sxy);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "SCALARS Syy float 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f\n", particles[i]->Syy);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "VECTORS velocity float\n");
	for (auto it = bodies.begin(); it != bodies.end(); ++it) {
		T **particles = it->get_particles();
		for (unsigned int i = 0; i < it->get_num_part(); i++) {
			fprintf(fp, "%f %f %f\n", particles[i]->vx, particles[i]->vy, 0.);
		}
	}
	fprintf(fp, "\n");

	fclose(fp);
}


#endif /* VTK_WRITER_H_ */
