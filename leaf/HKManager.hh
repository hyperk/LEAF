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

#define VERBOSE 		0

//WCSim parameters:
// If using a WCSim version without mPMT implementation
#define WCSIM_single_PMT_type
// If using a WCSim version with bugged number of PMT implementation
#define BUG_WCGEO

// Usefull processor functions:
#define GetPMTType(x)		x>=mPMT_ID_SHIFT?1:0
#define GetDistance(a,b)	sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) )
#define GetLength(a)		sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]) )
#define GetR(a)			sqrt( (a[0]*a[0]) + (a[1]*a[1]) )
#define GetScalarProd(a,b)	a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

// Hit definition:
#define NormalPMT		0 // B&L hit
#define MiniPMT			1 // 3" PMT hit
#define AllPMT			2

// Number of PMT configuration:
#define NPMT_CONFIGURATION 	2

// Hit parameters:
#define mPMT_ID_SHIFT 		 1000000
#define MAX_PMT 		10000000

// mPMT parameters:
#define mPMT_TOP_ID		19 //Number of 3" PMT per mPMT
//Different PMT config in an mPMT. I assumed here a rotational symetry of the mPMT, so PMT 1 to 12 are the same, 13 to 18 are the same and 19 is separated
#define NGROUP_PMT 		3

// Structure defining PMT hit:
struct PMTHit {
	double T;
	double Q;
	int PMT;
};

namespace HKGeoTools {

	void CrossProduct( double* a, double* b, double* c );
	void Normalize( double* a );
};


class HKManager/* : public TObject */{

	public:
		HKManager();
		~HKManager();
		static HKManager* GetME();
		
		void SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal=4200., double dDarkRate_mPMT=100. );
				
		// Manage hit info
		void ResetHitInfo() { fHitInfo.clear(); }
		void AddHit(double time, double charge, int pmtType, int tubeNumber) {
		
			// tubeNumber is from 1 to xxx in Hit array
			// but from 0 to xxx -1 in PMT info
			tubeNumber -= 1;
			
			if ( tubeNumber < 0 ) {
				std::cout << "ERROR: tubeNumber is below 0 (" << tubeNumber << ")" << std::endl;
			}
			struct PMTHit hHit;
	
			if ( pmtType == 0 ) { // Normal PMT
				hHit.PMT = tubeNumber;
			}
			else { 		      // mPMT
			
#ifdef WCSIM_single_PMT_type
				//Reject hit if it's not from a standard PMT
				return;
#endif

				hHit.PMT = mPMT_ID_SHIFT + tubeNumber;
			}
			
			hHit.T = time;
			hHit.Q = charge;
			
			fHitInfo.push_back(hHit);
		}
		
		int    GetNumberOfHit() 				{ return fHitInfo.size(); 		}
		PMTHit GetHitInfo(int iHit)				{ return fHitInfo[iHit]; 		}
		
		int    GetNPMT() 					{ return fNumberOfPMT; 		}
		double GetDarkNoise()					{ return fDarkRate_Normal;		}
		double GetDarkNoise_kHz()				{ return fDarkRate_Normal/1e3;	}
		
		int    GetNmPMT() 					{ return fNumberOfmPMT; 		}
		double GetDarkNoise_mPMT()				{ return fDarkRate_mPMT;		}
		double GetDarkNoise_mPMT_kHz()			{ return fDarkRate_mPMT/1e3;		}
		
		double GetMinimumDiagonal() 				{ return fPMTDiagDistance; 		}
		double GetBarrelExtremum_Lo()				{ return fPMTBarrelExtremZ_Lo;	}
		double GetBarrelExtremum_Hi()				{ return fPMTBarrelExtremZ_Hi;	}
		
		std::vector<double> GetPMTInfo(int iPMT) 		{ return fPMT_Info[iPMT]; 		}
		
		int GetMPMT_Reference(int iPMT) 			{ return fPMT_RefInMPMT[iPMT];	}
		std::vector<double> GetMPMT_RefX(int iPMT) 		{ return fPMT_RefX[iPMT];		}
		std::vector<double> GetMPMT_RefY(int iPMT) 		{ return fPMT_RefY[iPMT];		}
		std::vector<double> GetMPMT_RefZ(int iPMT) 		{ return fPMT_RefZ[iPMT];		}
		
		int GetMPMT_TubeID(int iPMT)				{ return fPMT_TubeInMPMT[iPMT]; 	}
		int GetMPMT_Group(int iPMT) 				{ return fPMT_Group[iPMT]; 		}
		
		double GetTankRadius()					{ return fTankRadius; 			}
		double GetTankHeight()					{ return fTankHeight; 			}
		double GetTankHalfHeight()				{ return fTankHalfHeight; 		}
				
	private:
	
		void LoadPMTInfo();
		void MakeMPMTReferencial(int iPMT);
		
		// Geometric function
		double GroupPMTs(int pmtType, int pmt_number);
		
	private:
		// WCSim objects:
		static HKManager* myManager;


		WCSimRootGeom * fWCGeo;
		double fDarkRate_Normal;
		double fDarkRate_mPMT;
		
		// Detector geometry:
		double fTankRadius;
		double fTankHeight;
		double fTankHalfHeight;
		
		// Darknoise
		double fDarkRate_ns[NPMT_CONFIGURATION];
		
		// PMT informations
		// Load every thing at the beginning to speed up fitting
		// mPMT tube ID start with 9000000
		std::vector< std::vector<double> > 	fPMT_Info; 
		// 0 -> PMT X
		// 1 -> PMT Y
		// 2 -> PMT Z
		// 3 -> PMT DIR X
		// 4 -> PMT DIR Y
		// 5 -> PMT DIR Z
		std::vector<int> 			fPMT_Group; 
		std::vector<int> 			fPMT_TubeInMPMT; 
		
		// mPMT Referencial
		std::vector<int> 			fPMT_RefInMPMT; 
		std::vector< std::vector<double> > 	fPMT_RefX;
		std::vector< std::vector<double> > 	fPMT_RefY;
		std::vector< std::vector<double> > 	fPMT_RefZ;
		
		// Information for reconstruction
		int 	fNumberOfPMT;
		int	fNumberOfmPMT;
		
		
		double 	fPMTDiagDistance;
		double 	fPMTBarrelExtremZ_Hi;
		double 	fPMTBarrelExtremZ_Lo;
		
		// Input/Output
		std::vector< PMTHit > fHitInfo;
};

#endif

