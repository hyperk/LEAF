#include "RootHitInfo.h"

RootHit::RootHit() {
	 
	PMT		= 0;
	PMT_original	= 0;
	
	T 		= 0.;
	Q 		= 0.;
	
}


RootHit::RootHit(int lID, TimeDelta::short_time_t lT, double lQ, bool mPMT) {

	PMT		= lID;
	PMT_original	= lID;
	
	if ( mPMT ) {
		PMT += HKAA::kmPMT_Shift;
	}
	
	T 		= lT;
	Q 		= lQ;

}

// Sort functor following Hit Time
bool RootHit::SortFunctor_HitTime::operator() (
		const RootHit &a,
		const RootHit &b) const {
		
	double ta = a.T;
	double tb = b.T;
	
	return ta < tb;
}


RootHitCollection::RootHitCollection() {
	timestamp = 0;
	first_unique = 0;
}

RootHitCollection::~RootHitCollection() {
	this->Clear();
}

void RootHitCollection::Clear() {
	timestamp = 0;
	first_unique = 0;
	hits.clear();	
}

void RootHitCollection::SortByTime() {
	std::sort(hits.begin(), hits.end(), RootHit::SortFunctor_HitTime());
}

bool RootHitCollection::Append(const RootHitCollection lHC) {
	// Adapted from SubSample append by T.Dealtry
	

	if ( hits.size() == 0 ) {
		// First time RootHitCollection is filled
		hits = lHC.hits;
		timestamp = lHC.timestamp;
	}
	else {
		// Need to shift the hit times by the difference of timestamp offsets
		TimeDelta::short_time_t time_shift = (lHC.timestamp - timestamp) / TimeDelta::ns;
	
		for ( unsigned int i = 0; i < lHC.Size(); i++ ) {
			
			RootHit lHit = lHC.hits[i];
			lHit.T += time_shift;
		
			hits.push_back( lHit );
		}
	}
	
	return true;
}
