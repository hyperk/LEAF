#ifndef HK_GEOMETRYPMT_V1_HPP
#define HK_GEOMETRYPMT_V1_HPP

#include "HKGeometryPMT.hpp"

/************************************************
 *  HKGeometryPMTV1 class definition
 *  Version 1 of HKGeometryPMT
 *
 * Just stores the PMT ID number used in software for this PMT.
 * Other information will be able to be looked up using channel mapping functionality (currently undefined).
 * This includes
 * * Position
 * * Orientation
 * * PMT type (20"/mPMT/OD PMT)
 * * PMT batch (e.g. 20" vX)
 * * PMT serial #
 * * Whether PMT is on barrel / top cap / bottom cap
 * * Electronics board type
 * * Electronics board serial #
 * * Electronics board channel #
 ***********************************************/
class HKGeometryPMTV1 : public HKGeometryPMT {
	public:

		HKGeometryPMTV1();  //!< Constructor

		//! Create HKGeometryPMTV1 from an HKGeometryPMT pointor
		/*!
		    \param in HKGeometryPMT*
		*/
		HKGeometryPMTV1(HKGeometryPMT* in);

		~HKGeometryPMTV1() override;  //!< Destructor
		void Reset() override;  //!< Reset function, called after each event. Use it to initialize variable

		static constexpr unsigned int GetVersionStatic() {
			return 1;
		}  //!< Static function. Return the class version number

		unsigned int GetVersion() const override {
			return HKGeometryPMTV1::GetVersionStatic();
		}  //!< Return the class version number

		// v1 functions

		//! Get PMTSoftwareID
		/*!
		    The PMT ID number used in software for this PMT
		    \return unsigned int
		*/
		unsigned int GetPMTSoftwareID() const override { return m_pmt_software_id; }

		//! Set PMTSoftwareID
		/*!
		    The PMT ID number used in software for this PMT
		    \param in unsigned int
		    \return `true` if successful, `false` if not
		*/
		bool SetPMTSoftwareID(unsigned int in) override {
			m_pmt_software_id = in;
			return true;
		}

		//! Get sub ID
		/*!
		    Return the PMT sub ID number. 0 for 20" PMT and OD PMT. 1-19 for mPMT's 3" PMTs
		    \return unsigned int
		*/
		unsigned int GetSubID() const override { return m_pmt_sub_id; }

		//! Set sub ID
		/*!
		    Set the PMT sub ID number. 0 for 20" PMT and OD PMT. 1-19 for mPMT's 3" PMTs
		    \param in unsigned int
		    \return `true` if successful, `false` if not
		*/
		bool SetSubID(unsigned int in) override {
			m_pmt_sub_id = in;
			return true;
		}

		//! Get PMT Type
		/*!
		    \return PMTType
		*/
		PMTType GetType() const override { return m_pmt_type; }

		//! Set PMT Type
		/*!
		    \param in PMTType
		    \return `true` if successful, `false` if not
		*/
		bool SetType(PMTType in) override {
			m_pmt_type = in;
			return true;
		}

		//! Get PMT Position in cm
		/*!
		    \return ROOT::Math::XYZVector. position over [x,y,z]
		*/
		virtual ROOT::Math::XYZVector GetPositionInCm() const override { return m_pmt_position; }

		//! Set PMT Position in cm
		/*!
		    \param x double
		    \param y double
		    \param z double
		    \return `true` if successful, `false` if not
		*/
		bool SetPositionInCm(double x, double y, double z) override {
			m_pmt_position.SetXYZ(x, y, z);
			return true;
		}

		//! Set PMT Position in cm
		/*!
		    \param x ROOT::Math::XYZVector
		    \return `true` if successful, `false` if not
		*/
		bool SetPositionInCm(ROOT::Math::XYZVector& in) override {
			m_pmt_position = in;
			return true;
		}

		//! Get PMT Orientation
		/*!
		    \return ROOT::Math::XYZVector. orientation over [x,y,z]
		*/
		ROOT::Math::XYZVector GetOrientation() const override { return m_pmt_orientation; }

		//! Set PMT Orientation
		/*!
		    \param x double
		    \param y double
		    \param z double
		    \return `true` if successful, `false` if not
		*/
		bool SetOrientation(double x, double y, double z) override {
			m_pmt_orientation.SetXYZ(x, y, z);
			return true;
		}

		//! Set PMT Orientation
		/*!
		    \param x ROOT::Math::XYZVector
		    \return `true` if successful, `false` if not
		*/
		bool SetPoSetOrientationsitionInCm(ROOT::Math::XYZVector& in) override {
			m_pmt_orientation = in;
			return true;
		}

	private:

		// v1 data members
		unsigned int m_pmt_software_id;        //!< PMTSoftwareID
		PMTType m_pmt_type;                    //!< PMT type
		unsigned int m_pmt_sub_id;             //!< PMT sub id (0 for 20" or OD PMT, 1-19 for mPMT's 3" PMTs)
		ROOT::Math::XYZVector m_pmt_position;  //!< Position in cm
		ROOT::Math::XYZVector m_pmt_orientation;  //!< Normalized orientation

		ClassDefOverride(HKGeometryPMTV1, 1);  //!< ROOT Class definition
};

#endif
