/*********************************************************************************/
/**	HKManager.hh								   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: June 17th 2020							   **/
/**	Desc: Manage fitter and PMT hits for HK analysis			   **/
/*********************************************************************************/


#ifndef HKManager_hh
#define HKManager_hh

//#include <algorithm>
//#include <ctime>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

//WCSim Headers
#include "WCSimRootGeom.hh"

//ROOT Headers
#include "TMath.h"

//DataModel informations
#include "HKGeometry.hpp"
#include "HKDarkNoise.hpp"
#include "HKGeometryPMT.hpp"
#include "HKHit.hpp"
#include "Constants/Suffix.hpp"
#include "Enums/PMTType.hpp"

//WCSim parameters:
// If using a WCSim version without mPMT implementation
//#define WCSIM_single_PMT_type

// Number of PMT configuration:
#define NPMT_CONFIGURATION 	1

class WCSimReader {

	public:
		WCSimReader();
		~WCSimReader();
	
#ifndef HK_USE_ROOT7	
		const HKGeometry* GetGeometry() { return fGeometry.get(); }
		const HKDarkNoise* GetDarkNoise() { return fDarkNoise.get(); }
		const HKGeometryPMTCollection* GetGeometryPMT_ID() { return fGeometryPMT_ID.get(); }
		const HKHitsCollection* GetHitCollection_ID() { return fHitCollection_ID.get(); }
		const HKGeometryPMTCollection* GetGeometryPMT_mPMT() { return fGeometryPMT_mPMT.get(); }
		const HKHitsCollection* GetHitCollection_mPMT() { return fHitCollection_mPMT.get(); }
#else
		std::shared_ptr<HKGeometry> GetGeometry() { return fGeometry; }
		std::shared_ptr<HKDarkNoise> GetDarkNoise() { return fDarkNoise; }
		std::shared_ptr<HKGeometryPMTCollection> GetGeometryPMT_ID() { return fGeometryPMT_ID; }
		std::shared_ptr<HKHitsCollection> GetHitCollection_ID() { return fHitCollection_ID; }
		std::shared_ptr<HKGeometryPMTCollection> GetGeometryPMT_mPMT() { return fGeometryPMT_mPMT; }
		std::shared_ptr<HKHitsCollection> GetHitCollection_mPMT() { return fHitCollection_mPMT; }
#endif
		void SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal=4200., double dDarkRate_mPMT=100. );
		
		// const HitCollection<Hit>* GetSecondaryHitCollection() { return &fSecondaryHitCollection; }
				
		// Manage hit info
		void ResetHitInfo() { fHitCollection_ID->Reset(); fHitCollection_mPMT->Reset(); }
		// void ResetSecondaryHitInfo() { fSecondaryHitCollection.Clean(); }

		void AddHit_ID(double time, double charge, int tubeNumber) {		
			// tubeNumber is from 1 to xxx in Hit array
			if ( tubeNumber < 1 ) {
				std::cout << "ERROR: tubeNumber is below 0 (" << tubeNumber << ")" << std::endl;
			}

			fHitCollection_ID->AddHit();
			fHitCollection_ID->back()->SetPMTNumber(tubeNumber);
			fHitCollection_ID->back()->SetTime(time);
			fHitCollection_ID->back()->SetCharge(charge);
		}

		void AddHit_mPMT(double time, double charge, int tubeNumber) {		
			// tubeNumber is from 1 to xxx in Hit array
			if ( tubeNumber < 1 ) {
				std::cout << "ERROR: tubeNumber is below 0 (" << tubeNumber << ")" << std::endl;
			}

			fHitCollection_mPMT->AddHit();
			fHitCollection_mPMT->back()->SetPMTNumber(tubeNumber);
			fHitCollection_mPMT->back()->SetTime(time);
			fHitCollection_mPMT->back()->SetCharge(charge);
		}

		// void AddSecondaryHit(double time, double charge, int pmtType, int tubeNumber) {
		
		// 	// tubeNumber is from 1 to xxx in Hit array
		// 	if ( tubeNumber < 1 ) {
		// 		std::cout << "ERROR: tubeNumber is below 0 (" << tubeNumber << ")" << std::endl;
		// 	}
			
		// 	Hit hHit (tubeNumber, time, charge, (HKAA::PMTType) pmtType);	
			
		// 	fSecondaryHitCollection.Add(hHit);
		// }
		
		/*
		// Commented on 2020/11/19 by Guillaume: I don't remember their purpose. To be removed?
		double GetMinimumDiagonal() 				{ return fPMTDiagDistance; 		}
		double GetBarrelExtremum_Lo()				{ return fPMTBarrelExtremZ_Lo;	}
		double GetBarrelExtremum_Hi()				{ return fPMTBarrelExtremZ_Hi;	}
		*/
			
	private:
		void LoadPMTInfo();
		
		// WCSim Geometry
		WCSimRootGeom * fWCGeo;

#ifndef HK_USE_ROOT7	
		// Geometry
		std::unique_ptr<HKGeometry> fGeometry;
		std::unique_ptr<HKDarkNoise> fDarkNoise;
		std::unique_ptr<HKGeometryPMTCollection> fGeometryPMT_ID;
		std::unique_ptr<HKGeometryPMTCollection> fGeometryPMT_mPMT;
		
		// HitCollection
		std::unique_ptr<HKHitsCollection> fHitCollection_ID;
		std::unique_ptr<HKHitsCollection> fHitCollection_mPMT;
#else
		// Geometry
		std::shared_ptr<HKGeometry> fGeometry;
		std::shared_ptr<HKDarkNoise> fDarkNoise;
		std::shared_ptr<HKGeometryPMTCollection> fGeometryPMT_ID;
		std::shared_ptr<HKGeometryPMTCollection> fGeometryPMT_mPMT;
		
		// HitCollection
		std::shared_ptr<HKHitsCollection> fHitCollection_ID;
		std::shared_ptr<HKHitsCollection> fHitCollection_mPMT;
#endif
				
		// PMT informations
		
		
		
		/*
		double 	fPMTDiagDistance;
		double 	fPMTBarrelExtremZ_Hi;
		double 	fPMTBarrelExtremZ_Lo;
		*/
		
};

class HKManager : public WCSimReader {

	public:
		HKManager();
		~HKManager();
		static HKManager* GetME();
		
	private:
		// WCSim objects:
		static HKManager* myManager;
	
};

#endif

