#ifndef RootTriggerInfo_H
#define RootTriggerInfo_H

#include "Environments.h"
#include "RootRecoInfo.h"
#include "TimeDelta.h"

#include "TObject.h"

/**
* \class RootTriggerInfo
 *
 * This class holds Trigger information 
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/16 14:20:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 

class RootTriggerInfo : public TObject {

	friend class EventInfo;

	public:
		RootTriggerInfo()			{				}
		
		unsigned int Id;
		HKAA::DAQType DAQ;
		
		TimeDelta time_start;
		TimeDelta time_end;
		TimeDelta time;
		
		int nhits;		
		
		// Reco Informations
		std::vector<RootRecoInfo*> RecoList;
		RootRecoInfo* RecoAt(int n) const	{ 	return RecoList[n];	}	
		
		
	private:
	
		RootTriggerInfo(const HKAA::DAQType lDAQ, unsigned int lId, TimeDelta lTime); 
		virtual ~RootTriggerInfo();
		
		ClassDef(RootTriggerInfo,1) //EventRootInfo structure
};

#if !defined(__CLING__)
ClassImp(RootTriggerInfo)
#endif

#endif
