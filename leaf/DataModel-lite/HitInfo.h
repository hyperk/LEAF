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

class Hit {


	public:
		Hit(); 
		Hit(int lID, TimeDelta::short_time_t lT, double lQ, bool mPMT=false); 

		int PMT;
		int PMT_original;
		TimeDelta::short_time_t T;
		double Q;
		
		
		// Sort functor following hit time after resolution applied
		struct SortFunctor_HitTime {
			bool operator() (const Hit &a, const Hit &b) const;
		};
		
};


class HitCollection {

	public:
		HitCollection(); 
		~HitCollection(); 
		
		int first_unique;
		TimeDelta timestamp;
		std::vector<Hit> hits;
		
		// Data access functions
		Hit &operator[](int n) {
			return hits[n];
		}
		
		// Information getter		
		unsigned int Size() const		{ return hits.size();		}
		Hit At(int n)	const			{ return hits[n];		}
		
		// Filler
		void Add(Hit lHit)			{ hits.push_back(lHit);	}
		bool Append(const HitCollection lHC);
		
		// Cleaner
		void Clear();
		
		void SortByTime();
};





#endif
