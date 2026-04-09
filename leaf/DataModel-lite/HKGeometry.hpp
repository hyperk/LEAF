#ifndef HK_GEOMETRY_HPP
#define HK_GEOMETRY_HPP

/*************************************************
 *  HKGeometry.hpp
 *   - Author:
 *   - Date: 1970/01/01
 ************************************************/

#include <array>
#include <string>

#include "Enums/DetectorID.hpp"
#include "HKGeometryPMT.hpp"
#include "HKObject.hpp"

/************************************************
 *  HKGeometry class definition
 ***********************************************/
class HKGeometry : public HKObject {
	public:

		HKGeometry();  //!< Constructor

		//! Create HKGeometry from an HKGeometry pointor
		/*!
		    \param in HKGeometry*
		*/
		HKGeometry(HKGeometry* in);

		virtual ~HKGeometry()                   = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the class version number

		static std::string GetBaseNameStatic() {
			return "HKGeometry";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKGeometry::GetBaseNameStatic();
		}  //!< Return the class base name

		// v1 functions
		//! Get IWCDElevatorHeight
		/*!
		    \return float
		*/
		virtual float GetIWCDElevatorHeight() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set IWCDElevatorHeight
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetIWCDElevatorHeight(UNUSED_PARAM float in) { return false; }

		//! Get IsCavernMaterialSimulated
		/*!
		    \return bool
		*/
		virtual bool GetIsCavernMaterialSimulated() const {
			return HKFormatError::ThrowErrorVersion<bool>(__PRETTY_FUNCTION__, 0);
		}

		//! Set IsCavernMaterialSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetIsCavernMaterialSimulated(UNUSED_PARAM bool in) { return false; }

		//! Get IsODSimulated
		/*!
		    \return bool
		*/
		virtual bool GetIsODSimulated() const {
			return HKFormatError::ThrowErrorVersion<bool>(__PRETTY_FUNCTION__, 0);
		}

		//! Set IsODSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetIsODSimulated(UNUSED_PARAM bool in) { return false; }

		//! Get IsDeadSpaceSimulated
		/*!
		    \return bool
		*/
		virtual bool GetIsDeadSpaceSimulated() const {
			return HKFormatError::ThrowErrorVersion<bool>(__PRETTY_FUNCTION__, 0);
		}

		//! Set IsDeadSpaceSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetIsDeadSpaceSimulated(UNUSED_PARAM bool in) { return false; }

		//! Get AreODPMTsActive
		/*!
		    \return bool
		*/
		virtual bool GetAreODPMTsActive() const {
			return HKFormatError::ThrowErrorVersion<bool>(__PRETTY_FUNCTION__, 0);
		}

		//! Set AreODPMTsActive
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetAreODPMTsActive(UNUSED_PARAM bool in) { return false; }

		//! Get OuterTyvekTotalHeight
		/*!
		    \return float
		*/
		virtual float GetOuterTyvekTotalHeight() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set OuterTyvekTotalHeight
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetOuterTyvekTotalHeight(UNUSED_PARAM float in) { return false; }

		//! Get OuterTyvekRadius
		/*!
		    \return float
		*/
		virtual float GetOuterTyvekRadius() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set OuterTyvekRadius
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetOuterTyvekRadius(UNUSED_PARAM float in) { return false; }

		//! Get InnerTyvekTotalHeight
		/*!
		    \return float
		*/
		virtual float GetInnerTyvekTotalHeight() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set InnerTyvekTotalHeight
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetInnerTyvekTotalHeight(UNUSED_PARAM float in) { return false; }

		//! Get InnerTyvekRadius
		/*!
		    \return float
		*/
		virtual float GetInnerTyvekRadius() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set InnerTyvekRadius
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetInnerTyvekRadius(UNUSED_PARAM float in) { return false; }

		//! Get BlacksheetTotalHeight
		/*!
		    \return float
		*/
		virtual float GetBlacksheetTotalHeight() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set BlacksheetTotalHeight
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetBlacksheetTotalHeight(UNUSED_PARAM float in) { return false; }

		//! Get BlacksheetRadius
		/*!
		    \return float
		*/
		virtual float GetBlacksheetRadius() const {
			return HKFormatError::ThrowErrorVersion<float>(__PRETTY_FUNCTION__, 0);
		}

		//! Set BlacksheetRadius
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetBlacksheetRadius(UNUSED_PARAM float in) { return false; }

		//! Get DetectorID
		/*!
		    \return DetectorID
		*/
		virtual DetectorID GetDetectorID() const {
			return HKFormatError::ThrowErrorVersion<DetectorID>(__PRETTY_FUNCTION__, DetectorID::kUndefined);
		}

		//! Set DetectorID
		/*!
		    \param in DetectorID
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetDetectorID(UNUSED_PARAM DetectorID in) { return false; }

		ClassDef(HKGeometry, 1);  //!< ROOT Class definition
};

#endif
