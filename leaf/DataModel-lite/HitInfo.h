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
		Hit(int lID, TimeDelta::short_time_t lT, double lQ, const HKAA::PMTType lPMT=HKAA::kIDPMT_BnL); 
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

#endif
