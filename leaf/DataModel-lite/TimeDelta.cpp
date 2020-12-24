#include <cmath>
#include "TimeDelta.h"

#ifdef ROOT5
const double TimeDelta::s_long_time_unit = 1.;
#endif

TimeDelta::TimeDelta(double naive_ns){
	m_long_time = 0;
	m_short_time = naive_ns;
	Normalize();
	// Avoid loss of precision when possible
	double difference = naive_ns - (*this / ns);
	m_short_time += difference;
	Normalize();
}

TimeDelta::~TimeDelta() {
}

void TimeDelta::Normalize(){
	long_time_t long_diff = std::floor(m_short_time / s_long_time_unit);
	m_long_time += long_diff;
	m_short_time -= (long_diff * s_long_time_unit);
}
/*
TimeDelta TimeDelta::operator=(const double& time) {
	TimeDelta new_delta(time);
	return new_delta;
}
*/
TimeDelta operator*(const TimeDelta& old_delta, double factor){
	TimeDelta new_delta(old_delta);
	new_delta.m_short_time *= factor;
	new_delta.m_long_time *= factor;
	// Make sure no information of m_long_time is lost
	if (factor != 0){
		double remainder = old_delta.m_long_time - new_delta.m_long_time * (1. / factor);
		new_delta.m_short_time += (remainder * TimeDelta::s_long_time_unit);
	}
	new_delta.Normalize();
	return new_delta;
}

TimeDelta operator*(double factor, const TimeDelta& old_delta){
	return old_delta * factor;
}

double operator/(const TimeDelta& left_delta, const TimeDelta& right_delta){
	double left_ns = (left_delta.m_long_time * (double)TimeDelta::s_long_time_unit);
	left_ns += left_delta.m_short_time;
	double right_ns = (right_delta.m_long_time * (double)TimeDelta::s_long_time_unit);
	right_ns += right_delta.m_short_time;
	return left_ns / right_ns;
}

TimeDelta operator+(const TimeDelta& left_delta, const TimeDelta& right_delta){
	TimeDelta new_delta(left_delta);
	new_delta.m_long_time += right_delta.m_long_time;
	new_delta.m_short_time += right_delta.m_short_time;
	new_delta.Normalize();
	return new_delta;
}

TimeDelta operator-(const TimeDelta& left_delta, const TimeDelta& right_delta){
	return left_delta + (-1.0 * right_delta);
}

TimeDelta& operator+=(TimeDelta& left_delta, const TimeDelta& right_delta){
	left_delta = left_delta + right_delta;
	return left_delta;
}

TimeDelta& operator-=(TimeDelta& left_delta, const TimeDelta& right_delta){
	left_delta = left_delta - right_delta;
	return left_delta;
}

bool operator==(const TimeDelta& left_delta, const TimeDelta& right_delta){
	TimeDelta A(left_delta);
	TimeDelta B(right_delta);
	A.Normalize();
	B.Normalize();
	return (A.m_long_time == B.m_long_time) and (A.m_short_time == B.m_short_time);
}

bool operator<(const TimeDelta& left_delta, const TimeDelta& right_delta){
	TimeDelta A(left_delta);
	TimeDelta B(right_delta);
	A.Normalize();
	B.Normalize();
	return (A.m_long_time < B.m_long_time)
			or ((A.m_long_time == B.m_long_time) and (A.m_short_time < B.m_short_time));
}

bool operator<=(const TimeDelta& left_delta, const TimeDelta& right_delta){
	return (left_delta < right_delta) or (left_delta == right_delta);
}

bool operator>(const TimeDelta& left_delta, const TimeDelta& right_delta){
	return not (left_delta <= right_delta);
}

bool operator>=(const TimeDelta& left_delta, const TimeDelta& right_delta){
	return not (left_delta < right_delta);
}

bool operator!=(const TimeDelta& left_delta, const TimeDelta& right_delta){
	return not (left_delta == right_delta);
}

std::ostream& operator<<(std::ostream& outs, const TimeDelta& delta){
	return outs << delta.m_short_time + delta.m_long_time*TimeDelta::s_long_time_unit << " ns";
}

// Unit constants
const TimeDelta TimeDelta::ps = TimeDelta(0.001);
const TimeDelta TimeDelta::ns = TimeDelta(1.0);
const TimeDelta TimeDelta::us = 1000 * TimeDelta::ns;
const TimeDelta TimeDelta::ms = 1000 * TimeDelta::us;
const TimeDelta TimeDelta::s = 1000 * TimeDelta::ms;
