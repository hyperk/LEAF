#ifndef HitCollection_H
#define HitCollection_H

#include "TimeDelta.h"
#include "HitInfo.h"

#include <algorithm>

/**
* \class HitCollection
 *
 * This class holds PMT hit information (PMT id, time, charge, etc.)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/12/09 16:15:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */

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
		void SetCollection(const HitCollection<T> lHC);
		void CopyCollection(const HitCollection<Hit> lHC);
		void CopyCollection(const HitCollection<HitExtended> lHC);
		void ResetTimes();
		
		void CopyDAQ(const HitCollection<T> lHC, const HKAA::DAQType iDAQ);
		void FilterTime(const HitCollection<T> lHC, const TimeDelta tStart, const TimeDelta tEnd);
		
		// Cleaner
		void Clean();
		void EraseFirstHit()			{ hits.erase(hits.begin());	}
		
		
		void SortByTime();
		void SortByTimeOfFlight();
};





#endif
