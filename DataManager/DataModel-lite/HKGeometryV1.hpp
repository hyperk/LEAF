#ifndef HK_GEOMETRY_V1_HPP
#define HK_GEOMETRY_V1_HPP

/*************************************************
 *  HKGeometryV1.hpp
 *   - Author:
 *   - Date: 1970/01/01
 ************************************************/

#include "HKGeometry.hpp"

/************************************************
 *  HKGeometryV1 class definition
 *  Version 1 of HKGeometry
 ***********************************************/
class HKGeometryV1 : public HKGeometry {
	public:

		HKGeometryV1();  //!< Constructor

		//! Create HKGeometryV1 from an HKGeometry pointor
		/*!
		    \param in HKGeometry*
		*/
		HKGeometryV1(HKGeometry* in);

		~HKGeometryV1() override;  //!< Destructor
		void Reset() override;     //!< Reset function, called after each event. Use it to initialize variable

		static constexpr unsigned int GetVersionStatic() {
			return 1;
		}  //!< Static function. Return the class version number

		unsigned int GetVersion() const override {
			return HKGeometryV1::GetVersionStatic();
		}  //!< Return the class version number

		// v1 functions
		//! Get IWCDElevatorHeight
		/*!
		    This is meaningless for detectors other than IWCD.
		    \return float
		*/
		float GetIWCDElevatorHeight() const override { return m_iwcd_elevator_height; }

		//! Set IWCDElevatorHeight
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetIWCDElevatorHeight(float in) override {
			m_iwcd_elevator_height = in;
			return true;
		}

		//! Get IsCavernMaterialSimulated
		/*!
		    \return bool
		*/
		bool GetIsCavernMaterialSimulated() const override { return m_is_cavern_material_simulated; }

		//! Set IsCavernMaterialSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		bool SetIsCavernMaterialSimulated(bool in) override {
			m_is_cavern_material_simulated = in;
			return true;
		}

		//! Get IsODSimulated
		/*!
		    \return bool
		*/
		bool GetIsODSimulated() const override { return m_is_od_simulated; }

		//! Set IsODSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		bool SetIsODSimulated(bool in) override {
			m_is_od_simulated = in;
			return true;
		}

		//! Get IsDeadSpaceSimulated
		/*!
		    \return bool
		*/
		bool GetIsDeadSpaceSimulated() const override { return m_is_dead_space_simulated; }

		//! Set IsDeadSpaceSimulated
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		bool SetIsDeadSpaceSimulated(bool in) override {
			m_is_dead_space_simulated = in;
			return true;
		}

		//! Get AreODPMTsActive
		/*!
		    \return bool
		*/
		bool GetAreODPMTsActive() const override { return m_are_od_pmts_active; }

		//! Set AreODPMTsActive
		/*!
		    \param in bool
		    \return `true` if successful, `false` if not
		*/
		bool SetAreODPMTsActive(bool in) override {
			m_are_od_pmts_active = in;
			return true;
		}

		//! Get the total height of the outer (cavern bottom/dome bottom side) tyvek
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetOuterTyvekTotalHeight() const override { return m_outer_tyvek_total_height; }

		//! Set the total height of the outer (cavern bottom/dome bottom side) tyvek
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetOuterTyvekTotalHeight(float in) override {
			m_outer_tyvek_total_height = in;
			return true;
		}

		//! Get the radius of the outer (cavern bottom/dome bottom side) tyvek
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetOuterTyvekRadius() const override { return m_outer_tyvek_radius; }

		//! Set the radius of the outer (cavern bottom/dome bottom side) tyvek
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetOuterTyvekRadius(float in) override {
			m_outer_tyvek_radius = in;
			return true;
		}

		//! Get the total height of the inner (PMT support structure side) tyvek
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetInnerTyvekTotalHeight() const override { return m_inner_tyvek_total_height; }

		//! Set the total height of the inner (PMT support structure side) tyvek
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetInnerTyvekTotalHeight(float in) override {
			m_inner_tyvek_total_height = in;
			return true;
		}

		//! Get the radius of the inner (PMT support structure side) tyvek
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetInnerTyvekRadius() const override { return m_inner_tyvek_radius; }

		//! Set the radius of the inner (PMT support structure side) tyvek
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetInnerTyvekRadius(float in) override {
			m_inner_tyvek_radius = in;
			return true;
		}

		//! Get the total height of the ID blacksheet
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetBlacksheetTotalHeight() const override { return m_blacksheet_total_height; }

		//! Set the total height of the ID blacksheet
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetBlacksheetTotalHeight(float in) override {
			m_blacksheet_total_height = in;
			return true;
		}

		//! Get the radius of the ID blacksheet
		/*!
		    Returns value in standard HK units: mm
		    \return float
		*/
		float GetBlacksheetRadius() const override { return m_blacksheet_radius; }

		//! Set the radius of the ID blacksheet
		/*!
		    \param in float
		    \return `true` if successful, `false` if not
		*/
		bool SetBlacksheetRadius(float in) override {
			m_blacksheet_radius = in;
			return true;
		}

		//! Get DetectorID
		/*!
		    \return DetectorID
		*/
		DetectorID GetDetectorID() const override { return m_detector_id; }

		//! Set DetectorID
		/*!
		    \param in DetectorID
		    \return `true` if successful, `false` if not
		*/
		bool SetDetectorID(DetectorID in) override {
			m_detector_id = in;
			return true;
		}

	private:

		// v1 data members
		float m_iwcd_elevator_height;         //!< IWCDElevatorHeight
		bool m_is_cavern_material_simulated;  //!< IsCavernMaterialSimulated
		bool m_is_od_simulated;               //!< IsODSimulated
		bool m_is_dead_space_simulated;       //!< IsDeadSpaceSimulated
		bool m_are_od_pmts_active;            //!< AreODPMTsActive
		float m_outer_tyvek_total_height;     //!< OuterTyvekTotalHeight
		float m_outer_tyvek_radius;           //!< OuterTyvekRadius
		float m_inner_tyvek_total_height;     //!< InnerTyvekTotalHeight
		float m_inner_tyvek_radius;           //!< InnerTyvekRadius
		float m_blacksheet_total_height;      //!< BlacksheetTotalHeight
		float m_blacksheet_radius;            //!< BlacksheetRadius
		DetectorID m_detector_id;             //!< DetectorID

		ClassDefOverride(HKGeometryV1, 1);  //!< ROOT Class definition
};

#endif
