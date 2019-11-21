/*
 * time_measurement.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: nguyenvu
 */
#include <boost/chrono.hpp>
#include <algorithm>
namespace bc =boost::chrono;
class time_measurement
{
	///start time
	bc::high_resolution_clock::time_point m_start;
	///stop time
	bc::high_resolution_clock::time_point m_stop;
	///result in microseconds
	double m_result = 0;
	public:
	void start()
	{
		m_start = bc::high_resolution_clock::now();
	}
	void stop()
	{
		m_stop = bc::high_resolution_clock::now();
		m_result = static_cast<double>(bc::duration_cast<bc::microseconds>(m_stop-m_start).count());
	}
	double measurement()
	{
		return m_result;
	}

};
