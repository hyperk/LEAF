#include "HKDarkNoiseV1.hpp"

HKDarkNoiseV1::HKDarkNoiseV1() {
	this->Reset();
}

HKDarkNoiseV1::HKDarkNoiseV1(HKDarkNoise* in) {
	// Load variable here:

	// v1 copies
	// m_a = in->GetA();
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
