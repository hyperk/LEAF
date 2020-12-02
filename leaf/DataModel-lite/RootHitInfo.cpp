#include "RootHitInfo.h"

RootHit::RootHit() {
	 
	PMT		= 0;
	PMT_original	= 0;
	
	T 		= 0.;
	Q 		= 0.;
	
	distance	= 0;
	ToF 		= 0;
		
	NormX		= 0;
	NormY		= 0;
	NormZ		= 0;
		
	theta		= 0;
	phi		= 0;
	
}


RootHit::~RootHit() {
}


RootHit::RootHit(int lID, TimeDelta::short_time_t lT, double lQ, bool mPMT) {

	PMT		= lID;
	PMT_original	= lID;
	
	if ( mPMT ) {
		PMT += HKAA::kmPMT_Shift;
	}
	
	T 		= lT;
	Q 		= lQ;
	
	distance	= 0;
	ToF 		= 0;
		
	NormX		= 0;
	NormY		= 0;
	NormZ		= 0;
		
	theta		= 0;
	phi		= 0;

}

// Sort functor following Hit Time
bool RootHit::SortFunctor_HitTime::operator() (
		const RootHit &a,
		const RootHit &b) const {
		
	double ta = a.T;
	double tb = b.T;
	
	return ta < tb;
}

// Sort functor following Hit Time
bool RootHit::SortFunctor_HitTimeOfFlight::operator() (
		const RootHit &a,
		const RootHit &b) const {
		
	double ta = a.ToF;
	double tb = b.ToF;
	
	return ta < tb;
}


RootHitCollection::RootHitCollection() {
	timestamp = 0;
	timestamp_last = 0;
	first_unique = 0;
}

RootHitCollection::~RootHitCollection() {
	this->Clear();
}

void RootHitCollection::Clean() {
	timestamp = 0;
	timestamp_last = 0;
	first_unique = 0;
	hits.clear();	
}

void RootHitCollection::SortByTime() {
	std::sort(hits.begin(), hits.end(), RootHit::SortFunctor_HitTime());
}
void RootHitCollection::SortByTimeOfFlight() {
	std::sort(hits.begin(), hits.end(), RootHit::SortFunctor_HitTimeOfFlight());
}

bool RootHitCollection::Append(const RootHitCollection lHC) {
	// Adapted from SubSample append by T.Dealtry

	if ( hits.size() == 0 ) {
		// First time RootHitCollection is filled
		hits = lHC.hits;
		timestamp = lHC.timestamp;
		timestamp_last = lHC.timestamp_last;
		first_unique = lHC.first_unique;
	}
	else {
		// Need to shift the hit times by the difference of timestamp offsets
		TimeDelta::short_time_t time_shift = (lHC.timestamp - timestamp) / TimeDelta::ns;
	
	
		for ( unsigned int i = 0; i < lHC.Size(); i++ ) {
			
			RootHit lHit = lHC.hits[i];
			lHit.T += time_shift;
		
			hits.push_back( lHit );
			
			timestamp_last = timestamp + TimeDelta(lHit.T);
		}
	}
	
	return true;
}

