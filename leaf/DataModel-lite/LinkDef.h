#ifdef __CINT__

	#pragma link off all globals;
	#pragma link off all classes;
	#pragma link off all functions;

	#pragma link C++ namespace HKFormatError + ;

	#pragma link C++ class vector < unsigned int> + ;
	#pragma link C++ class map < string, unsigned int> + ;

	#pragma link C++ class HKObject + ;
	#pragma link C++ class HKObjectCollection + ;

	#pragma link C++ class HKHit + ;
	#pragma link C++ class HKHitV1 + ;

	#pragma link C++ class HKHitsCollection + ;
	#pragma link C++ class vector < unique_ptr < HKHit, default_delete < HKHit>>> + ;
	#pragma link C++ class vector < shared_ptr < HKHit >> + ;
	#pragma link C++ class HKHitsCollectionT < HKHitV1> + ;

	#pragma link C++ class HKGeometryPMT + ;
	#pragma link C++ class HKGeometryPMTV1 + ;
	#pragma link C++ class HKGeometryPMTCollection + ;
	#pragma link C++ class vector < unique_ptr < HKGeometryPMT, default_delete < HKGeometryPMT>>> + ;
	#pragma link C++ class vector < shared_ptr < HKGeometryPMT >> + ;
	#pragma link C++ class HKGeometryPMTCollectionT < HKGeometryPMTV1> + ;

	#pragma link C++ class HKGeometry + ;
	#pragma link C++ class HKGeometryV1 + ;

	#pragma link C++ class HKDarkNoise+;
	#pragma link C++ class HKDarkNoiseV1+;

#endif
