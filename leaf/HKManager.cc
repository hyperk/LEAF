/*********************************************************************************/
/**	HKManager.cc							   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: June 17th 2020							   **/
/**	Desc: Manage fitter and PMT hits for HK analysis			   **/
/*********************************************************************************/

#include "HKManager.hh"

HKManager* HKManager::myManager=NULL;


/************************************************************************************************************************/
WCSimReader::WCSimReader() {
	fGeometry = new Geometry();
}

WCSimReader::~WCSimReader() {
	delete fGeometry;
}
/************************************************************************************************************************/

HKManager::HKManager() {
	myManager = this;
}

HKManager::~HKManager() {

	myManager = NULL;
	
	/*
	fPMTDiagDistance  = 0;
	fPMTBarrelExtremZ_Lo = 0;
	fPMTBarrelExtremZ_Hi = 0;
	*/
}

HKManager* HKManager::GetME() {

	if ( myManager ) return myManager;

	myManager = new HKManager();
	return myManager;
}

/************************************************************************************************************************/

void WCSimReader::SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal, double dDarkRate_mPMT ) {

	// Keep RootGeometry event:
	fWCGeo = wGeo;
	
	// Fill Geometry class
	fGeometry = new Geometry();
	
	fGeometry->detector_radius = fWCGeo->GetWCCylRadius();
	fGeometry->detector_length = fWCGeo->GetWCCylLength();
	
	fDarkRate_Normal = dDarkRate_Normal;
	fDarkRate_mPMT   = dDarkRate_mPMT;
	
	if ( fDarkRate_Normal == 0 ) {
#ifdef WCSIM_single_PMT_type	
		fDarkRate_ns[0] = 4200.  * fWCGeo->GetWCNumPMT() * 1e-9;
		fDarkRate_ns[1] = 0;
#else
		fDarkRate_ns[0] = 8400.  * fWCGeo->GetWCNumPMT(false) * 1e-9;
		fDarkRate_ns[1] = 100.   * fWCGeo->GetWCNumPMT(true ) * 1e-9 * HKAA::kmPMT_TopID;
#endif
	}
	else {	
#ifdef WCSIM_single_PMT_type	
		fDarkRate_ns[0] = fDarkRate_Normal * fWCGeo->GetWCNumPMT() * 1e-9;
		fDarkRate_ns[1] = 0;
#else

		int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
		int iNbr_mPMT = fWCGeo->GetWCNumPMT(true );

		fDarkRate_ns[0] = fDarkRate_Normal * iNbr_Norm * 1e-9;
		fDarkRate_ns[1] = fDarkRate_mPMT   * iNbr_mPMT * 1e-9 * HKAA::kmPMT_TopID;
#endif
	}
	
	fGeometry->pmt_dark_rate[HKAA::kIDPMT_BnL]	= fDarkRate_ns[0];
	fGeometry->pmt_dark_rate[HKAA::kIDPMT_3inch]	= fDarkRate_ns[1];
	

	this->LoadPMTInfo();
}

/************************************************************************************************************************/
// Geometry functions:

void WCSimReader::LoadPMTInfo() {

	// Get number of PMTs	
	
#ifdef WCSIM_single_PMT_type	
	int iNbr_Norm = fWCGeo->GetWCNumPMT();
	int iNbr_mPMT = 0;
#else
	int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
	int iNbr_mPMT = fWCGeo->GetWCNumPMT(true);
#endif

	// Fill geometry object	
	fGeometry->pmt_num[HKAA::kIDPMT_BnL]		= iNbr_Norm;
	fGeometry->pmt_num[HKAA::kIDPMT_3inch]	= iNbr_mPMT;


	std::cout << " LEAF setting # PMT = " << iNbr_Norm << std::endl;
	std::cout << " LEAF setting # mPMT = " << iNbr_mPMT << std::endl;
		
	WCSimRootPMT wPMT;
	/*
	double dMin1 = 0;
	double dMin2 = 0;
	*/
	
	// Normal PMTs
	for ( int iPMT=0; iPMT < iNbr_Norm; iPMT++ ) {
	
#ifdef WCSIM_single_PMT_type	
		wPMT = fWCGeo->GetPMT(iPMT);
#else
		wPMT = fWCGeo->GetPMT(iPMT,false);
#endif

		PMTInfo lPMTInfo;
		
		lPMTInfo.Id				= wPMT.GetTubeNo();
		
		lPMTInfo.Type				= HKAA::kIDPMT_BnL;
		lPMTInfo.mPMT				= false;
		
		for ( int j=0; j < 3; j++ ) {
			lPMTInfo.Position[j]		= wPMT.GetPosition(j);
			lPMTInfo.Orientation[j] 	= wPMT.GetOrientation(j);
		}
		/*
		if ( lPMTInfo.Position[2] < dMin1 ) {
			dMin2 = dMin1;
			dMin1 = lPMTInfo.Position[2];
		}
		
		if ( lPMTInfo.Position[2] > dMin1 && lPMTInfo.Position[2] < dMin2 ) {
			dMin2 = lPMTInfo.Position[2];
		}
		*/
		
		fGeometry->AddPMTInfo(lPMTInfo);
	}
	
	
	// Compute minimal diagonal distance:
	// PMT 4 is the closest diag PMT from PMT 0
	/*
	PMTInfo lInfoIdxZero = fGeometry->PMTList.at(HKAA::kID)[0];
	PMTInfo lInfoIdxFour = fGeometry->PMTList.at(HKAA::kID)[4];
	fPMTDiagDistance = std::ceil( sqrt( 	  pow(lInfoIdxZero.Position[0]-lInfoIdxFour.Position[0],2.) 
						+ pow(lInfoIdxZero.Position[1]-lInfoIdxFour.Position[1],2.) 
						+ pow(lInfoIdxZero.Position[2]-lInfoIdxFour.Position[2],2.) ) );
	fPMTDiagDistance += 2.; // Add 2cm for safety
	
	fPMTBarrelExtremZ_Lo = std::abs( std::ceil(dMin2) ) + 1;
	fPMTBarrelExtremZ_Hi = std::abs( std::ceil(dMin1) ) + 1;
	*/
	
#ifndef WCSIM_single_PMT_type	
	// mPMTs
	
	for ( int iPMT=0; iPMT < iNbr_mPMT; iPMT++ ) {
	
		wPMT = fWCGeo->GetPMT(iPMT,true);
		
		PMTInfo lPMTInfo;
		
		lPMTInfo.Id				= wPMT.GetTubeNo();
		lPMTInfo.Type				= HKAA::kIDPMT_3inch;
		lPMTInfo.mPMT				= true;
		
		for ( int j=0; j < 3; j++ ) {
			lPMTInfo.Position[j]		= wPMT.GetPosition(j);
			lPMTInfo.Orientation[j] 	= wPMT.GetOrientation(j);
		}
		
		lPMTInfo.mPMT_TubeNum			= wPMT.GetmPMT_PMTNo();
		 			
		if ( lPMTInfo.mPMT_TubeNum < 0 || lPMTInfo.mPMT_TubeNum > 30 ) {
			if ( iPMT > 0 ) {
				std::cout << " ERROR mPMT " << iPMT << "("<< lPMTInfo.Id << ") is not defined "
					<< "(" << lPMTInfo.Position[0] << ", " << lPMTInfo.Position[1] << ", " << lPMTInfo.Position[2] << ") " 
					<< "(" << lPMTInfo.Orientation[0] << ", " << lPMTInfo.Orientation[1] << ", " << lPMTInfo.Orientation[2] << ") " 
					<< lPMTInfo.mPMT_TubeNum << std::endl;
			}
			lPMTInfo.mPMT_TubeNum = 0;
			
		}
		fGeometry->AddPMTInfo(lPMTInfo);
	}
	
	fGeometry->Setup_mPMTs();
	
#endif
}
