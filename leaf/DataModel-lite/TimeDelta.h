#ifndef TIMEDELTA_H
#define TIMEDELTA_H

#include <iostream>
#include <vector>
#include <stdint.h>

/// Universal data type to store both large scale (unix timestamp) and short scale (hit times) time differences
///
/// The type and unit of the long and short scale times are optimised to allow
/// SubSamples to store many hits relative to a single timestamp efficiently.
///
/// The base unit for physics simulations (at least in WCSim) is ns, that that
/// is the chosen as the unit for the short time scales. A signed 64 bit
/// integer can hold values up to 2^63 = 9.22337204 x 10^18, which in ns would
/// be more than 290 years. Long enough for the long time deltas.
///
/// When creating TimeDelta instances, one should use the unit constants to do so:
///
///     TimeDelta t = 0.7 * TimeDelta::ms;
///
/// To convert into a naive floating point value, use them with the division operator:
///
///     double naive_t = t / TimeDelta::us;
///
class TimeDelta {

	public:
		/// Default constructor (sets all times to 0)
		TimeDelta() : m_short_time(0), m_long_time(0){};
		#ifndef ROOT5
		/// Copy constructor (just copies all member variables)
		TimeDelta(const TimeDelta&) = default;
		#else
		/// Copy constructor (just copies all member variables)
		TimeDelta(const TimeDelta& a) { m_short_time = a.m_short_time; m_long_time = a.m_long_time; }
		#endif
		/// Constructor from naive floating point value (in ns)
		TimeDelta(double naive_ns);
		
		~TimeDelta();

		/// Type for relative hit times within a SubSample. Unit = ns
		typedef double short_time_t;
		/// Type for absolute timestamps of a SubSample. Unit = ns
		typedef int64_t long_time_t;

		/// Member for short time deltas
		short_time_t m_short_time;
		/// Member for long time delta
		long_time_t m_long_time;
		
		#ifndef ROOT5
		/// Relative unit of long time member, i.e. long_unit / short_unit, both ns so = 1.
		static constexpr double s_long_time_unit = 1.;
		#else
		static const double s_long_time_unit;
		#endif

		/// Ensure that the time difference stored in m_short_time is small and positive.
		void Normalize();

		// Unit constants:

		/// TimeDelta of 1 ps
		static const TimeDelta ps;
		/// TimeDelta of 1 ns
		static const TimeDelta ns;
		/// TimeDelta of 1 us
		static const TimeDelta us;
		/// TimeDelta of 1 ms
		static const TimeDelta ms;
		/// TimeDelta of 1 s
		static const TimeDelta s;
		
		// Operators
		//TimeDelta operator=(const double& time);

};

	// Operators
	TimeDelta operator*(const TimeDelta& old_delta, double factor);
	TimeDelta operator*(double factor, const TimeDelta& old_delta);
	double operator/(const TimeDelta& left_delta, const TimeDelta& right_delta);
	TimeDelta operator+(const TimeDelta& left_delta, const TimeDelta& right_delta);
	TimeDelta operator-(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator<(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator<=(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator>(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator>=(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator==(const TimeDelta& left_delta, const TimeDelta& right_delta);
	bool operator!=(const TimeDelta& left_delta, const TimeDelta& right_delta);
	TimeDelta& operator+=(TimeDelta& left_delta, const TimeDelta& right_delta);
	TimeDelta& operator-=(TimeDelta& left_delta, const TimeDelta& right_delta);
	std::ostream& operator<<(std::ostream& outs, const TimeDelta& delta);


#endif
