#include "RootTriggerInfo.h"

RootTriggerInfo::RootTriggerInfo(unsigned int lId, TimeDelta lTime) {
	Id = lId;
	time = lTime;
	
	// Dummy value for now
	time_start	= time - 1e9;
	time_end	= time + 1e9;
		
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
