#ifndef SUFFIX_HPP
#define SUFFIX_HPP

#include <string>

// PMT collections

namespace HKSuffix {
#if HK_USE_ROOT7
	inline const std::string ID   = "ID";
	inline const std::string OD   = "OD";
	inline const std::string mPMT = "mPMT";
#else
	const std::string ID   = "ID";
	const std::string OD   = "OD";
	const std::string mPMT = "mPMT";
#endif

}  // namespace HKSuffix

#endif
