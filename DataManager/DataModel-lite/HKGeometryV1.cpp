#include "HKGeometryV1.hpp"

HKGeometryV1::HKGeometryV1() {
	this->Reset();
}

HKGeometryV1::HKGeometryV1(HKGeometry* in) {
	// Load variable here:

	// v1 copies
	m_iwcd_elevator_height         = in->GetIWCDElevatorHeight();
	m_is_cavern_material_simulated = in->GetIsCavernMaterialSimulated();
	m_is_od_simulated              = in->GetIsODSimulated();
	m_is_dead_space_simulated      = in->GetIsDeadSpaceSimulated();
	m_are_od_pmts_active           = in->GetAreODPMTsActive();
	m_outer_tyvek_total_height     = in->GetOuterTyvekTotalHeight();
	m_outer_tyvek_radius           = in->GetOuterTyvekRadius();
	m_inner_tyvek_total_height     = in->GetInnerTyvekTotalHeight();
	m_inner_tyvek_radius           = in->GetInnerTyvekRadius();
	m_blacksheet_total_height      = in->GetBlacksheetTotalHeight();
	m_blacksheet_radius            = in->GetBlacksheetRadius();
	m_detector_id                  = in->GetDetectorID();
}

HKGeometryV1::~HKGeometryV1() {
	this->Reset();
}

void HKGeometryV1::Reset() {
	// Initialize variable here:
	// Note, reset is called for every new event, so variables needing to be initialized
	// once for a full file should be defined in the creator

	// v1 resets
	m_iwcd_elevator_height         = 0;
	m_is_cavern_material_simulated = false;
	m_is_od_simulated              = false;
	m_is_dead_space_simulated      = false;
	m_are_od_pmts_active           = false;
	m_outer_tyvek_total_height     = 0;
	m_outer_tyvek_radius           = 0;
	m_inner_tyvek_total_height     = 0;
	m_inner_tyvek_radius           = 0;
	m_blacksheet_total_height      = 0;
	m_blacksheet_radius            = 0;
	m_detector_id                  = DetectorID::kUndefined;
}
