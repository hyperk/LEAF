#ifndef Environments_H
#define Environments_H

#include <string>
#include <stdint.h>
#include <vector>
#include <cmath>

/**
* \namespace Environments
 *
 * This namespace provides constants and enums for HK AstroAnalysis (lite version for LEAF)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/13 15:16:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *       
 */
 
#define mPMT_ID_SHIFT 		 1000000

 // Usefull processor functions:
#define GetPMTType(x)		x>=mPMT_ID_SHIFT?1:0
#define GetDistance(a,b)	sqrt( (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) )
#define GetLength(a)		sqrt( (a[0]*a[0]) + (a[1]*a[1]) + (a[2]*a[2]) )
#define GetR(a)		sqrt( (a[0]*a[0]) + (a[1]*a[1]) )
#define GetScalarProd(a,b)	a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

namespace HKAA {

	// Detector constants
	
	const int kMaxPMT_ID	= 10000000;
	const int kMaxPMT_OD	= 0;
	
	const int kmPMT_Shift	= 1000000;
	const int kmPMT_Groups	= 3;
	const int kmPMT_TopID	= 19; //Number of 3" PMT per mPMT
		
	// PMT Type
	enum PMTType {
		kIDPMT_BnL, 		// ID: Box and Line PMT
		kIDPMT_3inch, 		// ID: 3inch inside mPMT
		kODPMT_3inch,		// OD: 3inch 
		kPMTTypeMax
	};
	
	// DAQ Type
	enum DAQType {
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
}



#endif
