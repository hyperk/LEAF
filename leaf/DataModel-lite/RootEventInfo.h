#ifndef RootEventInfo_H
#define RootEventInfo_H

#include <vector>

#include "Environments.h"

#include "RootTriggerInfo.h"
#include "RootHitInfo.h"

#include "TObject.h"

/**
* \class EventInfo
 *
 * This class holds Event Information
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/16 13:42:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */

// Hold all triggers in one event
class RootEventInfo : public TObject {

	public:
		RootEventInfo();///< Simple constructor
		virtual ~RootEventInfo();

		std::vector<RootTriggerInfo*>	TriggerList;
		RootHitCollection		HitSample;
		
		HKAA::DAQType 			DAQLevel;

		ClassDef(RootEventInfo,1) //EventRootInfo structure
};

#if !defined(__CLING__)
ClassImp(RootEventInfo)
#endif

#endif
