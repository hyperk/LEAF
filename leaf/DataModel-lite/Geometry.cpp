#include "Geometry.h"
#include <iostream>

Geometry::Geometry(){

	// Initialize
	
	detector_length	= 0.;
	detector_radius	= 0.;
		
	HasOD			= false;
		
	pmt_radius		.resize(HKAA::kPMTTypeMax, 0.);
	pmt_num		.resize(HKAA::kPMTTypeMax, 0 );
	pmt_num_original	.resize(HKAA::kPMTTypeMax, 0 );
	pmt_dark_rate		.resize(HKAA::kPMTTypeMax, 0.);
	pmt_first		.resize(HKAA::kPMTTypeMax, HKAA::kMaxPMT_ID);
			
	PMTList[HKAA::kID].resize(HKAA::kMaxPMT_ID);
	PMTList[HKAA::kOD].resize(HKAA::kMaxPMT_OD);
}


Geometry::~Geometry(){

	PMTList[HKAA::kID].clear();
	PMTList[HKAA::kOD].clear();
}


void Geometry::AddPMTInfo(PMTInfo lInfo) {

	lInfo.Id_original	= lInfo.Id;
	lInfo.Masked		= false;

	if ( lInfo.Type == HKAA::kIDPMT_3inch ) {
		// Setup mPMT variables
		lInfo.Setup_mPMT();
	}
	
	if ( pmt_first[lInfo.Type] > lInfo.Id ) {
		pmt_first[lInfo.Type] = lInfo.Id;
	}
	
	if ( lInfo.Type == HKAA::kIDPMT_3inch || lInfo.Type == HKAA::kIDPMT_BnL ) {
		PMTList[HKAA::kID][lInfo.Id] = lInfo;
	}
	else if ( lInfo.Type == HKAA::kODPMT_3inch ) {
		PMTList[HKAA::kOD][lInfo.Id] = lInfo;
	}
	else {
		std::cout << "Unknown PMT Type " << lInfo.Type << std::endl;
	}
}

void Geometry::Setup_mPMTs() {

	for ( int i=HKAA::kmPMT_Shift; i < (HKAA::kmPMT_Shift + pmt_num[HKAA::kIDPMT_3inch]); i++ ) {
	
		PMTInfo lInfo	= PMTList[HKAA::kID][i];
		PMTInfo lInfoTop	= PMTList[HKAA::kID][lInfo.mPMT_RefTube];
		
		lInfo.Setup_mPMT_Referencial(lInfoTop);
	}
}
