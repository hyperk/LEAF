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
	pmt_first		.resize(HKAA::kPMTTypeMax, HKAA::kMaxPMT);
			
	PMTList.resize(HKAA::kMaxPMT);
}


Geometry::~Geometry(){

	PMTList.clear();
}


void Geometry::AddPMTInfo(PMTInfo lInfo) {

	lInfo.Id_original	= lInfo.Id;
	lInfo.Masked		= false;
	
	// Add Id shift:
	lInfo.Id += HKAA::kPMTId_Shift[lInfo.Type];;

	if ( lInfo.Type == HKAA::kIDPMT_3inch ) {
		// Setup mPMT variables
		lInfo.Setup_mPMT();
	}
	if ( pmt_first[lInfo.Type] > lInfo.Id ) {
		pmt_first[lInfo.Type] = lInfo.Id;
	}
	
	
	PMTList[lInfo.Id] = lInfo;
}

void Geometry::Finalise() {

	this->Setup_mPMTs();
	
	std::cout << "[DEBUG] PMT[" << HKAA::kIDPMT_BnL 	<< "] first Id is: " << pmt_first[HKAA::kIDPMT_BnL	] << ", Id shift is: " << HKAA::kPMTId_Shift[HKAA::kIDPMT_BnL  ] <<  std::endl;
	std::cout << "[DEBUG] PMT[" << HKAA::kIDPMT_3inch 	<< "] first Id is: " << pmt_first[HKAA::kIDPMT_3inch	] << ", Id shift is: " << HKAA::kPMTId_Shift[HKAA::kIDPMT_3inch] << std::endl;
	std::cout << "[DEBUG] PMT[" << HKAA::kODPMT_3inch 	<< "] first Id is: " << pmt_first[HKAA::kODPMT_3inch	] << ", Id shift is: " << HKAA::kPMTId_Shift[HKAA::kODPMT_3inch] << std::endl;
}

void Geometry::Setup_mPMTs() {

	for ( int i=HKAA::kmPMT_Shift; i < (HKAA::kmPMT_Shift + pmt_num[HKAA::kIDPMT_3inch]); i++ ) {
	
		PMTInfo lInfo		= PMTList[i];
		PMTInfo lInfoTop	= PMTList[lInfo.mPMT_RefTube];
		
		lInfo.Setup_mPMT_Referencial(lInfoTop);
	}
}
