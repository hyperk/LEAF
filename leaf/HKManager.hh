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
#include "Geometry.h"
#include "HitCollection.h"

//WCSim parameters:
// If using a WCSim version without mPMT implementation
//#define WCSIM_single_PMT_type

// Number of PMT configuration:
#define NPMT_CONFIGURATION 	2

class WCSimReader {

	public:
		WCSimReader();
		~WCSimReader();
		
		const Geometry* GetGeometry() { return fGeometry; }
		void SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal=4200., double dDarkRate_mPMT=100. );
		
		const HitCollection<Hit>* GetHitCollection() { return &fHitCollection; }
				
		// Manage hit info
		void ResetHitInfo() { fHitCollection.Clean(); }
		void AddHit(double time, double charge, int pmtType, int tubeNumber) {
		
			// tubeNumber is from 1 to xxx in Hit array
			if ( tubeNumber < 1 ) {
				std::cout << "ERROR: tubeNumber is below 0 (" << tubeNumber << ")" << std::endl;
			}
			
			Hit hHit (tubeNumber, time, charge, (HKAA::PMTType) pmtType);	
			
			fHitCollection.Add(hHit);
		}
		
		/*
		// Commented on 2020/11/19 by Guillaume: I don't remember their purpose. To be removed?
		double GetMinimumDiagonal() 				{ return fPMTDiagDistance; 		}
		double GetBarrelExtremum_Lo()				{ return fPMTBarrelExtremZ_Lo;	}
		double GetBarrelExtremum_Hi()				{ return fPMTBarrelExtremZ_Hi;	}
		*/
			
	private:
		void LoadPMTInfo();
		
		// Geometry
		WCSimRootGeom * fWCGeo;
		Geometry* fGeometry;
		
		// HitCollection
		HitCollection<Hit> fHitCollection;
		
		// Darknoise
		double fDarkRate_Normal;
		double fDarkRate_mPMT;
		
		// Darknoise
		double fDarkRate_ns[NPMT_CONFIGURATION];
		
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

