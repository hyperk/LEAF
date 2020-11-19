#ifndef Geometry_H
#define Geometry_H

#include <map>
#include <string>
#include <vector>

#include "Environments.h"

#include "PMTInfo.h"

/**
* \class Geometry
 *
 * This class holds Tank Geometry informations
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/16 09:37:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */

class Geometry {

	friend class WCSimReader;
	
	public:
		Geometry(); ///< Simple constructor
		
		// WCSim geometry information
		
		double detector_length;
		double detector_radius;
		
		bool HasOD;
		
		// PMT Info
		double	pmt_radius[HKAA::kPMTTypeMax];
		int 	pmt_num[HKAA::kPMTTypeMax];
		double	pmt_dark_rate[HKAA::kPMTTypeMax];
		
		std::map<const HKAA::DetectorType, std::vector<PMTInfo> > PMTList;
		
		const std::vector<PMTInfo> *GetPMTList(const HKAA::DetectorType lDetector=HKAA::kID) const { return &PMTList.at(lDetector); } ;
						
  		
	private: 
		// Function to properly set PMTInfo
		void AddPMTInfo(PMTInfo lInfo);
		
		// Start making referencial for mPMTs
		void Setup_mPMTs();
	 
};



#endif
