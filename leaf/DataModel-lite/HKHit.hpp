#ifndef HK_HIT_HPP
#define HK_HIT_HPP

/************************************************
 *  HKHit.hpp
 *   - Author: Guillaume Pronost
 *   - Date: 2024/07/10
 ***********************************************/

#include <memory>
#include <vector>

#include "HKObject.hpp"

/************************************************
 *  HKHit base class definition
 *
 *  Store PMT hit information
 ***********************************************/
class HKHit : public HKObject {
	public:

		HKHit();  //!< Constructor

		//! Create HKHit from an HKHit pointor
		/*!
		    \param in HKHit*
		*/
		HKHit(HKHit* in);

		virtual ~HKHit()                        = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the class version number

		static std::string GetBaseNameStatic() {
			return "HKHit";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKHit::GetBaseNameStatic();
		}  //!< Return the class base name

		// v1 functions
		//! Return PMT number
		/*!
		    \return unsigned short. PMT number
		*/
		virtual unsigned short GetPMTNumber() const {
			return HKFormatError::ThrowErrorVersion<unsigned short>(__PRETTY_FUNCTION__, 0);
		}

		//! Return hit's time
		//! TODO: Define unit
		/*!
		    \return double. Hit's time
		*/
		virtual double GetTime() const {
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, -9999.);
		}

		//! Return hit's charge
		//! TODO: Define unit
		/*!
		    \return double. Hit's charge
		*/
		virtual double GetCharge() const {
			return HKFormatError::ThrowErrorVersion<double>(__PRETTY_FUNCTION__, -9999.);
		}

		//! Set PMT number
		/*!
		    \param in unsigned short. PMT number
		*/
		virtual bool SetPMTNumber(UNUSED_PARAM unsigned short in) { return false; }

		//! Set hit's time
		//! TODO: Define unit
		/*!
		    \param in double. Hit's time
		*/
		virtual bool SetTime(UNUSED_PARAM double in) { return false; }

		//! Set hit's charge
		//! TODO: Define unit
		/*!
		    \param in double. Hit's charge
		*/
		virtual bool SetCharge(UNUSED_PARAM double in) { return false; }

		ClassDef(HKHit, 1);  //!< Class definition for ROOT
};

/************************************************
 *  HKHitsCollection base class definition
 *
 *  Collection class wrapping a vector of HKHit
 ***********************************************/
class HKHitsCollection : public HKObjectCollection {
	public:

		HKHitsCollection() {}  //!< Creator

		virtual ~HKHitsCollection()             = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the contained class version number

		static std::string GetBaseNameStatic() {
			return "HKHitsCollection";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKHitsCollection::GetBaseNameStatic();
		}  //!< Return the class base name

		//! operator []
		/*!
		    \param i integer. Entry index
		*/
		virtual HKHit* operator[](UNUSED_PARAM int i) {
			return HKFormatError::ThrowErrorVersion<HKHit*>(__PRETTY_FUNCTION__, NULL);
		}

		//! Wrapper for std::vector at()
		/*!
		    \param i integer. Entry index
		*/
		virtual HKHit* at(UNUSED_PARAM int i) {
			return HKFormatError::ThrowErrorVersion<HKHit*>(__PRETTY_FUNCTION__, NULL);
		}

		virtual const HKHit* at(UNUSED_PARAM int i) const {
			return HKFormatError::ThrowErrorVersion<HKHit*>(__PRETTY_FUNCTION__, NULL);
		}

		//! Wrapper for std::vector back()
		virtual HKHit* back() { return HKFormatError::ThrowErrorVersion<HKHit*>(__PRETTY_FUNCTION__, NULL); }

		//! Wrapper for std::vector size()
		virtual size_t size() const { return HKFormatError::ThrowErrorVersion<size_t>(__PRETTY_FUNCTION__, 0); }

		//! Wrapper for std::vector resize()
		/*!
		    \param entries integer. Number of entries
		*/
		virtual void resize(UNUSED_PARAM size_t entries) { return; }

		//! Wrapper for std::vector push_back()
		/*!
		    \param in HKHit*. Object to be inserted
		*/
		virtual void push_back(UNUSED_PARAM HKHit* in) { return; }

		//! Wrapper for std::vector emplace_back()
		/*!
		    \param in HKHit*. Object to be inserted
		*/
		virtual void emplace_back(UNUSED_PARAM HKHit* in) { return; }

		//! Emplace back an HKHit instance in the collection
		virtual void AddHit() { return; }

	private:

#ifndef HK_USE_ROOT7
		//! Return the pointor to the unique_ptr
		virtual std::unique_ptr<std::vector<std::unique_ptr<HKHit>>>* GetVectorUniquePointor() {
			return HKFormatError::ThrowErrorVersion<std::unique_ptr<std::vector<std::unique_ptr<HKHit>>>*>(
			    __PRETTY_FUNCTION__,
			    NULL);
		}

		//! Return the pointor to the vector
		virtual std::vector<std::unique_ptr<HKHit>>* GetVectorPointor() {
			return HKFormatError::ThrowErrorVersion<std::vector<std::unique_ptr<HKHit>>*>(__PRETTY_FUNCTION__,
			                                                                              NULL);
		}

		template<class> friend class HKBranchManagerCollection;
#endif

		ClassDef(HKHitsCollection, 1);  //!< Class definition for ROOT
};

/************************************************
 *  HKHitsCollectionTemplate definition
 *
 *  Collection class wrapping a vector of HKHit
 ***********************************************/
template<class T> class HKHitsCollectionT : public HKHitsCollection {

		static_assert(std::is_base_of<HKHit, T>::value, "T must inherit from HKHit");

	public:

		HKHitsCollectionT(bool initialize = true) {
			if(initialize) {
#ifdef HK_USE_ROOT7
				contents = std::make_shared<std::vector<std::unique_ptr<T>>>();
#else
				contents = std::make_unique<std::vector<std::unique_ptr<HKHit>>>();
#endif
				this->Reset();
			}
		}  //!< Creator

		//! Constructor. Initialize the collection with specified number of entries
		/*!
		    \param entries integer. Number of entries
		*/
		HKHitsCollectionT(size_t entries) {
			this->Reset();
			this->resize(entries);
		}

		~HKHitsCollectionT() { this->Reset(); }  //!< Destructor

		void Reset() override {
			contents->clear();
		}  //!< Reset function, called after each event. Clear the collection.

		static constexpr unsigned int GetVersionStatic() {
			return T::GetVersionStatic();
		}  //!< Static function. Return the contained class version number

		unsigned int GetVersion() const override {
			return T::GetVersionStatic();
		}  //!< Return the contained class version number

		//! operator []
		/*!
		    \param i integer. Entry index
		*/
		HKHit* operator[](int i) override { return (*contents)[i].get(); }  // Use at(i) here too?

		//! Wrapper for std::vector at()
		/*!
		    \param i integer. Entry index
		*/
		HKHit* at(int i) override { return contents->at(i).get(); }

		const HKHit* at(int i) const { return contents->at(i).get(); }

		//! Wrapper for std::vector back()
		HKHit* back() override { return contents->back().get(); }

		//! Wrapper for std::vector size()
		size_t size() const override { return contents->size(); }

		//! Wrapper for std::vector resize()
		/*!
		    \param entries integer. Number of entries
		*/
		void resize(size_t entries) override { 
#ifdef HK_USE_ROOT7
			for (auto& ptr : *contents) { ptr = std::make_unique<T>(); } 
#else
			for (auto& ptr : *contents) { ptr = std::make_unique<HKHit>(); } 
#endif
		}

		//! Wrapper for std::vector push_back()
		/*!
		    \param in HKHit*. Object to be inserted
		*/
		void push_back(HKHit* in) override { contents->push_back(std::make_unique<T>(in)); }

		//! Wrapper for std::vector emplace_back()
		/*!
		    \param in HKHit*. Object to be inserted
		*/
		void emplace_back(HKHit* in) override { contents->emplace_back(std::make_unique<T>(in)); }

		//! Emplace back an empty HKHit instance in the collection
		void AddHit() override { contents->emplace_back(std::make_unique<T>()); }

	private:

#ifdef HK_USE_ROOT7
		std::shared_ptr<std::vector<std::unique_ptr<T>>> contents;  //!< Vector of HKHit

		template<class, class, class> friend class HKDataMemberCollectionCreator;
#else
		//! Return the pointor to the unique_ptr
		std::unique_ptr<std::vector<std::unique_ptr<HKHit>>>* GetVectorUniquePointor() override {
			return &contents;
		}

		//! Return the pointor to the vector
		std::vector<std::unique_ptr<HKHit>>* GetVectorPointor() override { return contents.get(); }

		std::unique_ptr<std::vector<std::unique_ptr<HKHit>>> contents;  //!< Vector of HKHit

		template<class> friend class HKBranchManagerCollection;
#endif

		ClassDefOverride(HKHitsCollectionT, 1);  //!< Class definition for ROOT
};

#endif
