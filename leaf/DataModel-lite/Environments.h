#ifndef Environments_H
#define Environments_H

#include <string>
#include <stdint.h>
#include <vector>
#include <cmath>

/**
* \namespace Environments
 *
 * This namespace provides constants and enums for HK AstroAnalysis
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/13 15:16:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 
#define mPMT_ID_SHIFT 		 	100000
#define PMT_OD_SHIFT			 50000

 // Usefull processor functions:
#define Astro_GetPMTType(x)		x>=mPMT_ID_SHIFT?1:0

#define Astro_GetIfODPMT(x)		(x>=PMT_OD_SHIFT &&x<mPMT_ID_SHIFT)?1:0
#define Astro_GetIfmPMT(x)		(x>=mPMT_ID_SHIFT                 )?1:0
#define Astro_GetIfBnLPMT(x)		(                  x<PMT_OD_SHIFT )?1:0

#define Astro_GetDistance(a,b)	sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) )
#define Astro_GetLength(a)		sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]) )
#define Astro_GetR(a)			sqrt( (a[0]*a[0]) + (a[1]*a[1]) )
#define Astro_GetScalarProd(a,b)	a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

namespace HKAA {

	// Vertex reconstruction algorithm
	enum RecoVtxInfo {
		kRecoVtxNone = -1, 	//
		kRecoVtxMax
	};
	
	// Energy reconstruction algorithm
	enum RecoEneInfo {
		kRecoEneNone = -1, 	//
		kRecoEneMax
	};
	
	// Direction reconstruction algorithm
	enum RecoDirInfo {
		kRecoDirNone = -1, 	//
		kRecoDirMax
	};
	
	// Detector constants

	
	// PMT Type
	enum PMTType {
		kIDPMT_BnL, 		// ID: Box and Line PMT
		kIDPMT_3inch, 		// ID: 3inch inside mPMT
		kODPMT_3inch,		// OD: 3inch 
		kPMTTypeMax
	};
	
	const int kPMTId_Shift[kPMTTypeMax] = {
		0, 			// ID - Box and Line PMTs
		mPMT_ID_SHIFT, 	// ID - mPMTs
		PMT_OD_SHIFT 		// OD - PMTs
	};
	
	// DAQ Type
	enum DAQType {
		kDAQNone = -1, 	// 
		kIDDAQ_BnL, 		// ID: Box And Line PMT
		kIDDAQ_3inch, 		// ID: multi-PMT
		kIDDAQ_Hybrid,		// ID: merge B&L and mPMT 
		kODDAQ,		// OD 
		kDAQTypeMax
	};
	
	// Detector Type
	enum DetectorType {
		kID,			// ID 
		kOD,			// OD 
		kDetectorTypeMax
	};
	
	// Trigger Type
	// From WCSim Enumeration
	enum TriggerType {
		kTriggerUndefined = -1,	// Default 
		kTriggerNDigits,
		kTriggerNDigitsTest,
		kTriggerTestVertices,
		kTriggerNoTrig,
		kTriggerFailure,
		kTriggerTypeMax
	};	
	
	const int kMaxPMT	= 1000000;
	
	const int kmPMT_Shift	= kPMTId_Shift[kIDPMT_3inch];
	const int kmPMT_Groups	= 3;
	const int kmPMT_TopID	= 19; //Number of 3" PMT per mPMT
	
	const int kODPMT_Shift	= kPMTId_Shift[kODPMT_3inch];
	
}

#endif
