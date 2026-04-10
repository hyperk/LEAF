/*********************************************************************************/
/**	HKManager.cc							   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: June 17th 2020							   **/
/**	Desc: Manage fitter and PMT hits for HK analysis			   **/
/*********************************************************************************/

#include "HKManager.hpp"

HKManager* HKManager::myManager=NULL;

#include "HKGeometryV1.hpp"
#include "HKGeometryPMTV1.hpp"
#include "HKDarkNoiseV1.hpp"
#include "HKHitV1.hpp"

/************************************************************************************************************************/
WCSimReader::WCSimReader() {
#ifndef HK_USE_ROOT7
	fGeometry = std::make_unique<HKGeometryV1>();
	fDarkNoise = std::make_unique<HKDarkNoiseV1>();
	fGeometryPMT_ID = std::make_unique<HKGeometryPMTCollectionT<HKGeometryPMTV1>>();
	fGeometryPMT_mPMT = std::make_unique<HKGeometryPMTCollectionT<HKGeometryPMTV1>>();

	fHitCollection_ID = std::make_unique<HKHitsCollectionT<HKHitV1>>();
	fHitCollection_mPMT = std::make_unique<HKHitsCollectionT<HKHitV1>>();
#else
	fGeometry = std::static_pointer_cast<HKGeometry>(std::make_shared<HKGeometryV1>());
	fDarkNoise = std::static_pointer_cast<HKDarkNoise>(std::make_shared<HKDarkNoiseV1>());

	fGeometryPMT_ID = std::static_pointer_cast<HKGeometryPMTCollection>(std::make_shared<HKGeometryPMTCollectionT<HKGeometryPMTV1>>());
	fGeometryPMT_mPMT = std::static_pointer_cast<HKGeometryPMTCollection>(std::make_shared<HKGeometryPMTCollectionT<HKGeometryPMTV1>>());
	
	fHitCollection_ID = std::static_pointer_cast<HKHitsCollection>(std::make_shared<HKHitsCollectionT<HKHitV1>>());
	fHitCollection_mPMT = std::static_pointer_cast<HKHitsCollection>(std::make_shared<HKHitsCollectionT<HKHitV1>>());
#endif
}

WCSimReader::~WCSimReader() {
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
		
	fGeometry->SetInnerTyvekRadius(fWCGeo->GetWCCylRadius());
	fGeometry->SetInnerTyvekTotalHeight(fWCGeo->GetWCCylLength());
		
	if ( dDarkRate_Normal == 0 ) {
#ifdef WCSIM_single_PMT_type	
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kID, 	4200.  * fWCGeo->GetWCNumPMT() * 1e-9);
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kmPMT, 0.0);
#else
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kID, 	8400.  * fWCGeo->GetWCNumPMT(false) * 1e-9);
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kmPMT, 100.   * fWCGeo->GetWCNumPMT(true ) * 1e-9 * 19); // 19 is the number of 3" PMTs per mPMT
#endif
	}
	else {	
#ifdef WCSIM_single_PMT_type	
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kID,   dDarkRate_Normal * fWCGeo->GetWCNumPMT() * 1e-9);
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kmPMT, 0.0);
#else

		int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
		int iNbr_mPMT = fWCGeo->GetWCNumPMT(true );

		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kID,   dDarkRate_Normal * iNbr_Norm * 1e-9);
		fDarkNoise->SetAverageTotalDarkNoisePerNS(PMTType::kmPMT, dDarkRate_mPMT   * iNbr_mPMT * 1e-9 * 19); // 19 is the number of 3" PMTs per mPMT
#endif
	}	

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
	fGeometryPMT_ID->resize(iNbr_Norm ? iNbr_Norm+1 : 0);
	fGeometryPMT_mPMT->resize(iNbr_mPMT? iNbr_mPMT+1 : 0);

	std::cout << " LEAF setting # PMT = " << iNbr_Norm << std::endl;
	std::cout << " LEAF setting # mPMT = " << iNbr_mPMT << std::endl;
		
	WCSimRootPMT wPMT;
	/*
	double dMin1 = 0;
	double dMin2 = 0;
	*/
	
	// Normal PMTs
	std::cout << "Filling ID PMTs" << std::endl;
	for ( int iPMT=0; iPMT < iNbr_Norm; iPMT++ ) {	
#ifdef WCSIM_single_PMT_type	
		wPMT = fWCGeo->GetPMT(iPMT);
#else
		wPMT = fWCGeo->GetPMT(iPMT,false);
#endif

		int tubeNo = wPMT.GetTubeNo();
		/*
		std::cout << iPMT << " " 
					<< "TubeNo: " << wPMT.GetTubeNo() << " "
					<< "Pos: (" << wPMT.GetPosition(0) << ", " << wPMT.GetPosition(1) << ", " << wPMT.GetPosition(2) << ") "
					<< "Dir: (" << wPMT.GetOrientation(0) << ", " << wPMT.GetOrientation(1) << ", " << wPMT.GetOrientation(2) << ") "
					<< std::endl;

		std::cout << " check " << fGeometryPMT_ID->at(iPMT)->GetPMTSoftwareID() << std::endl;
		*/
		fGeometryPMT_ID->at(tubeNo)->SetPMTSoftwareID(wPMT.GetTubeNo());
		fGeometryPMT_ID->at(tubeNo)->SetType(PMTType::kID);
		fGeometryPMT_ID->at(tubeNo)->SetPositionInCm(wPMT.GetPosition(0), wPMT.GetPosition(1), wPMT.GetPosition(2));
		fGeometryPMT_ID->at(tubeNo)->SetOrientation(wPMT.GetOrientation(0), wPMT.GetOrientation(1), wPMT.GetOrientation(2));

		/*
		if ( wPMT.GetPosition(2) < dMin1 ) {
			dMin2 = dMin1;
			dMin1 = wPMT.GetPosition(2);
		}
		
		if ( wPMT.GetPosition(2) > dMin1 && wPMT.GetPosition(2) < dMin2 ) {
			dMin2 = wPMT.GetPosition(2);
		}
		*/
	}	
	
	// Compute minimal diagonal distance:
	// PMT 4 is the closest diag PMT from PMT 0
	/*
	PMTInfo lInfoIdxZero = fGeometry->PMTList.at(PMTType::kID)[0];
	PMTInfo lInfoIdxFour = fGeometry->PMTList.at(PMTType::kID)[4];
	fPMTDiagDistance = std::ceil( sqrt( 	  pow(lInfoIdxZero.Position[0]-lInfoIdxFour.Position[0],2.) 
						+ pow(lInfoIdxZero.Position[1]-lInfoIdxFour.Position[1],2.) 
						+ pow(lInfoIdxZero.Position[2]-lInfoIdxFour.Position[2],2.) ) );
	fPMTDiagDistance += 2.; // Add 2cm for safety
	
	fPMTBarrelExtremZ_Lo = std::abs( std::ceil(dMin2) ) + 1;
	fPMTBarrelExtremZ_Hi = std::abs( std::ceil(dMin1) ) + 1;
	*/
	
#ifndef WCSIM_single_PMT_type	
	// mPMTs	
	std::cout << "Filling mPMT" << std::endl;
	for ( int iPMT=0; iPMT < iNbr_mPMT; iPMT++ ) {	
		wPMT = fWCGeo->GetPMT(iPMT,true);
		int tubeNo = wPMT.GetTubeNo();

		fGeometryPMT_mPMT->at(tubeNo)->SetPMTSoftwareID(wPMT.GetTubeNo());
		fGeometryPMT_mPMT->at(tubeNo)->SetType(PMTType::kmPMT);
		fGeometryPMT_mPMT->at(tubeNo)->SetPositionInCm(wPMT.GetPosition(0), wPMT.GetPosition(1), wPMT.GetPosition(2));
		fGeometryPMT_mPMT->at(tubeNo)->SetOrientation(wPMT.GetOrientation(0), wPMT.GetOrientation(1), wPMT.GetOrientation(2));
		fGeometryPMT_mPMT->at(tubeNo)->SetSubID(wPMT.GetmPMT_PMTNo());
	}	
#endif
}
