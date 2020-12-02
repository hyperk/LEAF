#ifndef HitInfo_H
#define HitInfo_H

#include "Environments.h"
#include "TimeDelta.h"

#include <algorithm>

/**
* \class HitInfo
 *
 * This class holds PMT hit information (PMT id, time, charge, etc.)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/13 16:15:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */

class RootHit : public TObject {


	public:
		RootHit(); 
		RootHit(int lID, TimeDelta::short_time_t lT, double lQ, bool mPMT=false); 
		~RootHit(); 

		int PMT;
		int PMT_original;
		TimeDelta::short_time_t T;
		double Q;
		
		// Need for AstroAnalysis
		double distance;
		double ToF;
		
		double NormX;
		double NormY;
		double NormZ;
		
		double theta;
		double phi;
				
		// Sort functor following hit time after resolution applied
		struct SortFunctor_HitTime {
			bool operator() (const RootHit &a, const RootHit &b) const;
		};
		// Sort functor following hit time of flight after resolution applied
		struct SortFunctor_HitTimeOfFlight {
			bool operator() (const RootHit &a, const RootHit &b) const;
		};
		
		ClassDef(RootHit,1) //EventRootInfo structure
		
};


class RootHitCollection : public TObject {

	public:
		RootHitCollection(); 
		virtual ~RootHitCollection(); 
		
		int first_unique;
		TimeDelta timestamp;
		TimeDelta timestamp_last;
		std::vector<RootHit> hits;
		
		// Data access functions
		RootHit &operator[](int n) 		{ return hits[n];		}
		
		// Information getter		
		unsigned int Size() const		{ return hits.size();		}
		RootHit At(int n) const		{ return hits[n];		}
		
		// Filler
		void Add(RootHit lHit)			{ hits.push_back(lHit);	}
		bool Append(const RootHitCollection lHC);
		
		// Cleaner
		void Clean();
		void EraseFirstHit()			{ hits.erase(hits.begin());	}
		
		void SortByTime();
		void SortByTimeOfFlight();
				
		ClassDef(RootHitCollection,1) //EventRootInfo structure
};

#if !defined(__CLING__)
ClassImp(RootHitCollection)
ClassImp(RootHit)
#endif




#endif
