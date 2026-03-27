#ifndef HK_DARKNOISE_V1_HPP
#define HK_DARKNOISE_V1_HPP

/*************************************************
 *  HKDarkNoiseV1.hpp
 *   - Author:
 *   - Date: 1970/01/01
 ************************************************/

#include "HKDarkNoise.hpp"

/************************************************
 *  HKDarkNoiseV1 class definition
 *  Version 1 of HKDarkNoise
 ***********************************************/
class HKDarkNoiseV1 : public HKDarkNoise {
	public:

		HKDarkNoiseV1();  //!< Constructor

		//! Create HKDarkNoiseV1 from an HKDarkNoise pointor
		/*!
		    \param in HKDarkNoise*
		*/
		HKDarkNoiseV1(HKDarkNoise* in);

		~HKDarkNoiseV1() override;  //!< Destructor
		void Reset() override;    //!< Reset function, called after each event. Use it to initialize variable

		static constexpr unsigned int GetVersionStatic() {
			return 1;
		}  //!< Static function. Return the class version number

		unsigned int GetVersion() const override {
			return HKDarkNoiseV1::GetVersionStatic();
		}  //!< Return the class version number

		// v1 functions

		//! Return average total number of dark noise hit in 1ns for a given PMT type 
		//! Note: the value is ~ N_PMT * darknoise_per_PMT * 1e-9
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		*/
		double GetAverageTotalDarkNoisePerNS(PMTType in) const override { 
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, 0.0); 
		}
		
		//! Set average total number of dark noise hit in 1ns for a given PMT type 
		//! Note: the value is ~ N_PMT * darknoise_per_PMT * 1e-9
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param noise double. Average total dark noise hit
		*/
		bool SetAverageTotalDarkNoisePerNS(PMTType in, double noise) override {
			return false;
		}

		//! Return dark noise value for a given PMT in Hz
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param num int. Channel number 
		*/
		double GetDarkNoiseHz(PMTType in, int num) const override { 
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, 0.0); 
		}
		
		//! Set dark noise value for a given PMT type in Hz
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param num int. Channel number 
		    \param noise double. Dark noise value in Hz
		*/
		bool SetDarkNoiseHz(PMTType in, int num, double noise) override {
			return false;
		}

		// define function here:
		// int GetA() const override 		{ return m_a; }
		// bool SetA(int in) override		{ m_a = in; return true; }

	private:

		// v1 data members
		std::array<double, static_cast<size_t>(PMTType::kNumPMTTypes)> m_avg_dark_noise;  // Average dark noise for each PMT type (ID, OD, mPMT)
		std::array<std::vector<double>, static_cast<size_t>(PMTType::kNumPMTTypes)> m_dark_noise;  // Dark noise for each PMT (ID, OD, mPMT). 

		ClassDefOverride(HKDarkNoiseV1, 1);  //!< ROOT Class definition
};

#endif
