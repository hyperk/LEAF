#include "HKHitV1.hpp"

HKHitV1::HKHitV1() {
	this->Reset();
}

HKHitV1::HKHitV1(HKHit* in) {
	m_pmt_number = in->GetPMTNumber();
	m_time       = in->GetTime();
	m_charge     = in->GetCharge();
}

HKHitV1::~HKHitV1() {
	this->Reset();
}

void HKHitV1::Reset() {
	m_pmt_number = 0;
	m_time       = 0;
	m_charge     = 0;
}
