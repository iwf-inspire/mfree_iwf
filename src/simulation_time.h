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

//singleton to hold simulation time

#ifndef SIMULATION_TIME_H_
#define SIMULATION_TIME_H_

template <typename T> class rk4;

class simulation_time {
	template <typename T> friend class rk4;
public:
	static simulation_time& getInstance();
	simulation_time(simulation_time const &) = delete;
	void operator=(simulation_time const &) = delete;
	double get_time() const;
	double get_dt() const;

private:
	simulation_time();
	double m_time;
	double m_dt;
	void increment_time();
};

#endif /* SIMULATION_TIME_H_ */
