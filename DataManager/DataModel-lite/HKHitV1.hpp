#ifndef HK_HIT_V1_HPP
#define HK_HIT_V1_HPP

/************************************************
 *  HKHitV1.hpp
 *   - Author: Guillaume Pronost
 *   - Date: 2024/07/10
 ***********************************************/

#include "HKHit.hpp"

/************************************************
 *  HKHitV1 class definition
 *
 *  Version 1 of HKHit
 *
 *  Store PMT hit information
 ***********************************************/
class HKHitV1 : public HKHit {
	public:

		HKHitV1();  //!< Creator

		//! Create HKHitV1 from an HKHit pointor
		/*!
		    \param in HKHit*
		*/
		HKHitV1(HKHit* in);

		~HKHitV1() override;    //!< Destructor
		void Reset() override;  //!< Reset function, called after each event

		static constexpr unsigned int GetVersionStatic() {
			return 1;
		}  //!< Static function. Return the class version number

		unsigned int GetVersion() const override {
			return HKHitV1::GetVersionStatic();
		}  //!< Return the class version number

		// v1 functions
		//! Return PMT number
		/*!
		    \return unsigned short. PMT number
		*/
		unsigned short GetPMTNumber() const override { return m_pmt_number; }

		//! Return hit's time
		//! TODO: Define unit
		/*!
		    \return double. Hit's time
		*/
		double GetTime() const override { return m_time; }

		//! Return hit's charge
		//! TODO: Define unit
		/*!
		    \return double. Hit's charge
		*/
		double GetCharge() const override { return m_charge; }

		//! Set PMT number
		/*!
		    \param in unsigned short. PMT number
		*/
		bool SetPMTNumber(unsigned short in) override {
			m_pmt_number = in;
			return true;
		}

		//! Set hit's time
		//! TODO: Define unit
		/*!
		    \param in double. Hit's time
		*/
		bool SetTime(double in) override {
			m_time = in;
			return true;
		}

		//! Set hit's charge
		//! TODO: Define unit
		/*!
		    \param in double. Hit's charge
		*/
		bool SetCharge(double in) override {
			m_charge = in;
			return true;
		}

	private:

		unsigned short m_pmt_number;  //!< PMT number
		double m_time;                //!< Hit's time
		double m_charge;              //!< Hit's charge

		ClassDefOverride(HKHitV1, 1);  //!< Class definition for ROOT
};

#endif
