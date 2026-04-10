#include "HKDarkNoiseV1.hpp"

HKDarkNoiseV1::HKDarkNoiseV1() {
	m_avg_dark_noise[static_cast<size_t>(PMTType::kID)]   = 0.0;
	m_avg_dark_noise[static_cast<size_t>(PMTType::kOD)]   = 0.0;
	m_avg_dark_noise[static_cast<size_t>(PMTType::kmPMT)] = 0.0;

	this->Reset();
}

HKDarkNoiseV1::HKDarkNoiseV1(HKDarkNoise* in) {
	// Load variable here:

	// v1 copies
	m_avg_dark_noise[static_cast<size_t>(PMTType::kID)]   = in->GetAverageTotalDarkNoisePerNS(PMTType::kID);
	m_avg_dark_noise[static_cast<size_t>(PMTType::kOD)]   = in->GetAverageTotalDarkNoisePerNS(PMTType::kOD);
	m_avg_dark_noise[static_cast<size_t>(PMTType::kmPMT)] = in->GetAverageTotalDarkNoisePerNS(PMTType::kmPMT);

	m_dark_noise[static_cast<size_t>(PMTType::kID)]   = in->GetDarkNoiseVectorHz(PMTType::kID);
	m_dark_noise[static_cast<size_t>(PMTType::kOD)]   = in->GetDarkNoiseVectorHz(PMTType::kOD);
	m_dark_noise[static_cast<size_t>(PMTType::kmPMT)] = in->GetDarkNoiseVectorHz(PMTType::kmPMT);
}

HKDarkNoiseV1::~HKDarkNoiseV1() {
	this->Reset();
}

void HKDarkNoiseV1::Reset() {
	// Initialize variable here:
	// Note, reset is called for every new event, so variables needing to be initialized
	// once for a full file should be defined in the creator

	// v1 resets
	// m_a = 0;
}
