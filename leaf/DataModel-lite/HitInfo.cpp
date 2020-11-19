#include "HitInfo.h"

Hit::Hit() {
	 
	PMT		= 0;
	PMT_original	= 0;
	
	T 		= 0.;
	Q 		= 0.;
	
}


Hit::Hit(int lID, TimeDelta::short_time_t lT, double lQ, bool mPMT) {

	PMT		= lID;
	PMT_original	= lID;
	
	if ( mPMT ) {
		PMT += HKAA::kmPMT_Shift;
	}
	
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


HitCollection::HitCollection() {
	timestamp = 0;
	first_unique = 0;
}

HitCollection::~HitCollection() {
	this->Clear();
}

void HitCollection::Clear() {
	timestamp = 0;
	first_unique = 0;
	hits.clear();	
}

void HitCollection::SortByTime() {
	std::sort(hits.begin(), hits.end(), Hit::SortFunctor_HitTime());
}

bool HitCollection::Append(const HitCollection lHC) {
	// Adapted from SubSample append by T.Dealtry
	

	if ( hits.size() == 0 ) {
		// First time HitCollection is filled
		hits = lHC.hits;
		timestamp = lHC.timestamp;
	}
	else {
		// Need to shift the hit times by the difference of timestamp offsets
		TimeDelta::short_time_t time_shift = (lHC.timestamp - timestamp) / TimeDelta::ns;
	
		for ( unsigned int i = 0; i < lHC.Size(); i++ ) {
			
			Hit lHit = lHC.hits[i];
			lHit.T += time_shift;
		
			hits.push_back( lHit );
		}
	}
	
	return true;
}
