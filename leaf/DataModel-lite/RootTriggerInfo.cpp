#include "RootTriggerInfo.h"

RootTriggerInfo::RootTriggerInfo(const HKAA::DAQType lDAQ, unsigned int lId, TimeDelta lTime) {
	Id = lId;
	time = lTime;
	DAQ = lDAQ;
	
	// Dummy value for now
	time_start	= time;
	time_end	= time;
		
	RecoList.resize(HKAA::kRecoVtxMax);
	for ( unsigned int i = 0; i < RecoList.size(); i++ ) {
		RecoList[i] = new RootRecoInfo();
	}
	
}

RootTriggerInfo::~RootTriggerInfo() {
	for ( unsigned int i = 0; i < RecoList.size(); i++ ) {
		delete RecoList[i];
	}
	RecoList.clear();
}
