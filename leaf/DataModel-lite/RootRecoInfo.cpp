#include "RootRecoInfo.h"


RootRecoInfo::RootRecoInfo()
{	
	AlgVtx = HKAA::kRecoVtxNone; 
	
	Reset();
}

RootRecoInfo::~RootRecoInfo()
{	
	// Clear
	Vtx.clear();
	E.clear();
	Dir.clear();
	Dir_Goodness.clear();
}


void RootRecoInfo::Reset()
{
	// Clear
	Vtx.clear();
	E.clear();
	Dir.clear();
	Dir_Goodness.clear();
	
	ExtendedHitCollection.Clean();
	InTime20.Clean();
	InTime30.Clean();
	InTime50.Clean();
	
	// Re-initialize
	std::vector<double> lEmptyThree(3,0.);
	std::vector<double> lEmptyTwo(2,0.);
	
	Goodness 	= 0;	
	Goodness_Time 	= 0;
	
	Wall		= 0;
	
	Vtx		.resize(4,0.);		
	E		.resize(HKAA::kRecoEneMax,0.);
	Dir		.resize(HKAA::kRecoDirMax,lEmptyThree);
	Cone		.resize(HKAA::kRecoDirMax,lEmptyTwo);
	Dir_Goodness	.resize(HKAA::kRecoDirMax,0);
}


