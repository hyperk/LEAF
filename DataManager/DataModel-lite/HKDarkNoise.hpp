#ifndef HK_DARKNOISE_HPP
#define HK_DARKNOISE_HPP

/*************************************************
 *  HKDarkNoise.hpp
 *   - Author:
 *   - Date: 1970/01/01
 ************************************************/

#include <array>
#include <string>

#include "Enums/PMTType.hpp"
#include "HKObject.hpp"

/************************************************
 *  HKDarkNoise class definition
 ***********************************************/
class HKDarkNoise : public HKObject {
	public:

		HKDarkNoise();  //!< Constructor

		//! Create HKDarkNoise from an HKDarkNoise pointor
		/*!
		    \param in HKDarkNoise*
		*/
		HKDarkNoise(HKDarkNoise* in);

		virtual ~HKDarkNoise()                  = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the class version number

		static std::string GetBaseNameStatic() {
			return "HKDarkNoise";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKDarkNoise::GetBaseNameStatic();
		}  //!< Return the class base name

		// v1 functions

		//! Initialize dark noise vector for a given PMT type
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param number int. Number of PMTs
		*/
		virtual bool Initialize(PMTType in, int number) { return false; }

		//! Return average total number of dark noise hit in 1ns for a given PMT type
		//! Note: the value is ~ N_PMT * darknoise_per_PMT * 1e-9
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		*/
		virtual double GetAverageTotalDarkNoisePerNS(PMTType in) const {
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, 0.0);
		}

		//! Set average total number of dark noise hit in 1ns for a given PMT type
		//! Note: the value is ~ N_PMT * darknoise_per_PMT * 1e-9
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param noise double. Average total dark noise hit
		*/
		virtual bool SetAverageTotalDarkNoisePerNS(PMTType in, double noise) { return false; }

		//! Return dark noise rate for a given PMT in Hz
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param num int. Channel number
		*/
		virtual double GetDarkNoiseHz(PMTType in, int num) const {
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, 0.0);
		}

		//! Set dark noise value for a given PMT type in Hz
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param num int. Channel number
		    \param noise double. Dark noise value in Hz
		*/
		virtual bool SetDarkNoiseHz(PMTType in, int num, double noise) { return false; }

		//! Return dark noise vector for a given PMT type
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \return std::vector<double>. Dark noise values in Hz
		*/
		virtual const std::vector<double> GetDarkNoiseVectorHz(PMTType in) const { 
			return HKFormatError::ThrowErrorVersion<std::vector<double>>(__PRETTY_FUNCTION__, std::vector<double>());
		}

		//! Set dark noise vector for a given PMT type
		/*!
		    \param in PMTType. PMT type (ID, OD, mPMT)
		    \param noise std::vector<double>. Dark noise values in Hz
		*/
		virtual bool SetDarkNoiseVectorHz(PMTType in, const std::vector<double>& noise) { return false; }

		ClassDef(HKDarkNoise, 1);  //!< ROOT Class definition
};

#endif
