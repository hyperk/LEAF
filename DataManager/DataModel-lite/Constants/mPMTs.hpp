#ifndef MPMTS_HPP
#define MPMTS_HPP

/************************************************
 *  Constant related to mPMTs enum class definition
 ***********************************************/

namespace mPMTConstants {
#if HK_USE_ROOT7
	inline const int kNGroups = 3;
	inline const int kNPMTs   = 19; //Number of 3" PMT per mPMT
#else
	const int kNGroups = 3;
	const int kNPMTs   = 19; //Number of 3" PMT per mPMT
#endif

}  // namespace mPMTConstants

#endif