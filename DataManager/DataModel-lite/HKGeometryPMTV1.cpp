#include "HKGeometryPMTV1.hpp"

HKGeometryPMTV1::HKGeometryPMTV1() {

	m_pmt_software_id = 0;
	m_pmt_sub_id      = 0;
	m_pmt_position.SetXYZ(-9999., -9999., -9999.);
	m_pmt_orientation.SetXYZ(-9999., -9999., -9999.);
	m_pmt_bad_flag = false;

	this->Reset();
}

HKGeometryPMTV1::HKGeometryPMTV1(HKGeometryPMT* in) {
	// Load variable here:

	// v1 copies
	m_pmt_software_id = in->GetPMTSoftwareID();
	m_pmt_sub_id      = in->GetSubID();
	m_pmt_position    = in->GetPositionInCm();
	m_pmt_orientation = in->GetOrientation();
	m_pmt_bad_flag    = in->GetBadFlag();
}

HKGeometryPMTV1::~HKGeometryPMTV1() {
	this->Reset();
}

void HKGeometryPMTV1::Reset() {
	// Initialize variable here:
	// Note, reset is called for every new event, so variables needing to be initialized
	// once for a full file should be defined in the creator

	// v1 resets
}
