#include "HitInfo.h"


/*************************************************************************/

Hit::Hit() {
	 
	PMT		= 0;
	PMT_original	= 0;
	
	T 		= 0.;
	Q 		= 0.;
}


Hit::~Hit() { }


Hit::Hit(int lID, TimeDelta::short_time_t lT, double lQ, const HKAA::PMTType lPMT) {

	PMT		= lID + HKAA::kPMTId_Shift[lPMT];
	PMT_original	= lID;
	
	T 		= lT;
	Q 		= lQ;
}

// Sort functor following Hit Time
bool Hit::SortFunctor_HitTime::operator() (
		const Hit &a,
		const Hit &b) const {
		
	double ta = a.T;
	double tb = b.T;
	
	return ta < tb;
}

/*************************************************************************/

HitExtended::HitExtended() {

	distance	= 0;
	ToF 		= 0;
		
	NormX		= 0;
	NormY		= 0;
	NormZ		= 0;
		
	theta		= 0;
	phi		= 0;
}

HitExtended::HitExtended(const Hit& tHit) {


	PMT 		= tHit.PMT;
	PMT_original	= tHit.PMT_original;
	T		= tHit.T;
	Q		= tHit.Q;

	distance	= 0;
	ToF 		= 0;
		
	NormX		= 0;
	NormY		= 0;
	NormZ		= 0;
		
	theta		= 0;
	phi		= 0;
}

// Sort functor following Hit Time
bool HitExtended::SortFunctor_HitTimeOfFlight::operator() (
		const HitExtended &a,
		const HitExtended &b) const {
		
	double ta = a.ToF;
	double tb = b.ToF;
	
	return ta < tb;
}

