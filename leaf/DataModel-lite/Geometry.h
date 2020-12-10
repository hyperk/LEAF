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
	friend class PMTMask;
	
	public:
		Geometry(); ///< Simple constructor
		~Geometry();
		
		// WCSim geometry information
		
		double detector_length;
		double detector_radius;
		
		bool HasOD;
		
		// PMT Info
		std::vector<double>	pmt_radius;
		std::vector<int> 	pmt_num;
		std::vector<int> 	pmt_num_original;
		std::vector<double>	pmt_dark_rate;
		std::vector<int> 	pmt_first;
		
		std::map<const HKAA::DetectorType, std::vector<PMTInfo> > PMTList;
		
		const std::vector<PMTInfo> *GetPMTList(const HKAA::DetectorType lDetector=HKAA::kID) const { return &PMTList.at(lDetector); } ;
		
		PMTInfo GetPMT(int lChannel, const HKAA::DetectorType lDetector=HKAA::kID)	const {	return PMTList.at(lDetector)[lChannel];		}
		bool GetIfMasked(int lChannel, const HKAA::DetectorType lDetector=HKAA::kID)		const { 	return PMTList.at(lDetector)[lChannel].Masked;	}
  		
	private: 
		// Function to properly set PMTInfo
		void AddPMTInfo(PMTInfo lInfo);
		
		// Start making referencial for mPMTs
		void Setup_mPMTs();
		
		// Mask PMT
		void Mask_PMT(int lChannel, const HKAA::DetectorType lDetector=HKAA::kID)	{	PMTList.at(lDetector)[lChannel].Masked = true;		}
		void UnMask_PMT(int lChannel, const HKAA::DetectorType lDetector=HKAA::kID)	{	PMTList.at(lDetector)[lChannel].Masked = false;		}
};

#endif
