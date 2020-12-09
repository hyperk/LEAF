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
		~Hit(); 

		int PMT;
		int PMT_original;
		TimeDelta::short_time_t T;
		double Q;
						
		// Sort functor following hit time after resolution applied
		struct SortFunctor_HitTime {
			bool operator() (const Hit &a, const Hit &b) const;
		};
		
};

class HitExtended : public Hit {


	public:
		HitExtended();
		HitExtended(const Hit&);
		
		double distance;
		double ToF;
		
		double NormX;
		double NormY;
		double NormZ;
		
		double theta;
		double phi;
				
		// Sort functor following hit time of flight after resolution applied
		struct SortFunctor_HitTimeOfFlight {
			bool operator() (const HitExtended &a, const HitExtended &b) const;
		};
		
};


template <typename T> class HitCollection {

	public:
		HitCollection(); 
		virtual ~HitCollection(); 
		
		int first_unique;
		TimeDelta timestamp;
		TimeDelta timestamp_last;
		std::vector<T> hits;
		
		// Data access functions
		T &operator[](int n) 			{ return hits[n];		}
		
		// Information getter		
		unsigned int Size() const		{ return hits.size();		}
		T At(int n) const			{ return hits[n];		}
		
		// Filler
		void Add(T lHit)			{ hits.push_back(lHit);	}
		void Append(const HitCollection<T> lHC);
		void CopyCollection(const HitCollection<Hit> lHC);
		
		// Cleaner
		void Clean();
		void EraseFirstHit()			{ hits.erase(hits.begin());	}
		
		
		void SortByTime();
		void SortByTimeOfFlight();
};





#endif
