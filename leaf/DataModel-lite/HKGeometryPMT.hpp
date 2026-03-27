#ifndef HK_GEOMETRYPMT_HPP
#define HK_GEOMETRYPMT_HPP

/*************************************************
 *  HKGeometryPMT.hpp
 *   - Author: T. Dealtry
 *   - Date: 2024/11/14
 ************************************************/

#include <array>
#include <memory>
#include <string>

#include "Enums/PMTType.hpp"
#include "HKObject.hpp"

/************************************************
 *  HKGeometryPMT class definition
 ***********************************************/
class HKGeometryPMT : public HKObject {
	public:

		HKGeometryPMT();  //!< Constructor

		//! Create HKGeometryPMT from an HKGeometryPMT pointor
		/*!
		    \param in HKGeometryPMT*
		*/
		HKGeometryPMT(HKGeometryPMT* in);

		virtual ~HKGeometryPMT()                = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the class version number

		static std::string GetBaseNameStatic() {
			return "HKGeometryPMT";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKGeometryPMT::GetBaseNameStatic();
		}  //!< Return the class base name

		// v1 functions
		//! Get PMTSoftwareID
		/*!
		    \return unsigned int
		*/
		virtual unsigned int GetPMTSoftwareID() const {
			return HKFormatError::ThrowErrorVersion<unsigned int>(__PRETTY_FUNCTION__, 0);
		}

		//! Set PMTSoftwareID
		/*!
		    \param in unsigned int
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetPMTSoftwareID(UNUSED_PARAM unsigned int in) { return false; }

		//! Get PMT Type
		/*!
		    \return PMTType
		*/
		virtual PMTType GetType() const {
			return HKFormatError::ThrowErrorVersion<PMTType>(__PRETTY_FUNCTION__, PMTType::kUndefined);
		}

		//! Set PMT Type
		/*!
		    \param in PMTType
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetType(UNUSED_PARAM PMTType in) { return false; }

		//! Get sub ID
		/*!
		    Return the PMT sub ID number. 0 for 20" PMT and OD PMT. 1-19 for mPMT's 3" PMTs
		    \return unsigned int
		*/
		virtual unsigned int GetSubID() const {
			return HKFormatError::ThrowErrorVersion<unsigned int>(__PRETTY_FUNCTION__, 0);
		}

		//! Set sub ID
		/*!
		    Set the PMT sub ID number. 0 for 20" PMT and OD PMT. 1-19 for mPMT's 3" PMTs
		    \param in unsigned int
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetSubID(UNUSED_PARAM unsigned int in) { return false; }

		//! Get PMT Position in cm
		/*!
		    \return ROOT::Math::XYZVector. position over [x,y,z]
		*/
		virtual ROOT::Math::XYZVector GetPositionInCm() const {
			return HKFormatError::ThrowErrorVersion<ROOT::Math::XYZVector>(__PRETTY_FUNCTION__,
			                                                               ROOT::Math::XYZVector());
		}

		//! Set PMT Position in cm
		/*!
		    \param x double
		    \param y double
		    \param z double
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetPositionInCm(UNUSED_PARAM double x, UNUSED_PARAM double y, UNUSED_PARAM double z) {
			return false;
		}

		//! Set PMT Position in cm
		/*!
		    \param x ROOT::Math::XYZVector
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetPositionInCm(UNUSED_PARAM ROOT::Math::XYZVector& in) { return false; }

		//! Get PMT Orientation
		/*!
		    \return ROOT::Math::XYZVector. orientation over [x,y,z]
		*/
		virtual ROOT::Math::XYZVector GetOrientation() const {
			return HKFormatError::ThrowErrorVersion<ROOT::Math::XYZVector>(__PRETTY_FUNCTION__,
			                                                               ROOT::Math::XYZVector());
		}

		//! Set PMT Orientation
		/*!
		    \param x double
		    \param y double
		    \param z double
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetOrientation(UNUSED_PARAM double x, UNUSED_PARAM double y, UNUSED_PARAM double z) {
			return false;
		}

		//! Set PMT Orientation
		/*!
		    \param x ROOT::Math::XYZVector
		    \return `true` if successful, `false` if not
		*/
		virtual bool SetPoSetOrientationsitionInCm(UNUSED_PARAM ROOT::Math::XYZVector& in) { return false; }

		ClassDef(HKGeometryPMT, 1);  //!< ROOT Class definition
};

/************************************************
 *  HKGeometryPMTCollection base class definition
 *
 *  Collection class wrapping a vector of HKGeometryPMT
 ***********************************************/
class HKGeometryPMTCollection : public HKObjectCollection {
	public:

		HKGeometryPMTCollection() {}  //!< Creator

		virtual ~HKGeometryPMTCollection()      = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the contained class version number

		static std::string GetBaseNameStatic() {
			return "HKGeometryPMTCollection";
		}  //!< Static function. Return the class base name

		std::string GetBaseName() const {
			return HKGeometryPMTCollection::GetBaseNameStatic();
		}  //!< Return the class base name

		//! operator []
		/*!
		    \param i integer. Entry index
		*/
		virtual HKGeometryPMT* operator[](UNUSED_PARAM int i) {
			return HKFormatError::ThrowErrorVersion<HKGeometryPMT*>(__PRETTY_FUNCTION__, NULL);
		}

		//! Wrapper for std::vector at()
		/*!
		    \param i integer. Entry index
		*/
		virtual HKGeometryPMT* at(UNUSED_PARAM int i) {
			return HKFormatError::ThrowErrorVersion<HKGeometryPMT*>(__PRETTY_FUNCTION__, NULL);
		}

		virtual const HKGeometryPMT* at(UNUSED_PARAM int i) const {
			return HKFormatError::ThrowErrorVersion<HKGeometryPMT*>(__PRETTY_FUNCTION__, NULL);
		}


		//! Wrapper for std::vector back()
		virtual HKGeometryPMT* back() {
			return HKFormatError::ThrowErrorVersion<HKGeometryPMT*>(__PRETTY_FUNCTION__, NULL);
		}

		//! Wrapper for std::vector size()
		virtual size_t size() const { return HKFormatError::ThrowErrorVersion<size_t>(__PRETTY_FUNCTION__, 0); }

		//! Wrapper for std::vector resize()
		/*!
		    \param entries integer. Number of entries
		*/
		virtual void resize(UNUSED_PARAM size_t entries) { return; }

		//! Wrapper for std::vector push_back()
		/*!
		    \param in HKGeometryPMT*. Object to be inserted
		*/
		virtual void push_back(UNUSED_PARAM HKGeometryPMT* in) { return; }

		//! Wrapper for std::vector emplace_back()
		/*!
		    \param in HKGeometryPMT*. Object to be inserted
		*/
		virtual void emplace_back(UNUSED_PARAM HKGeometryPMT* in) { return; }

		//! Emplace back an HKGeometryPMT instance in the collection
		virtual void AddGeometryPMT() { return; }

	private:

#ifndef HK_USE_ROOT7
		//! Return the pointor to the unique_ptr
		virtual std::unique_ptr<std::vector<std::unique_ptr<HKGeometryPMT> > >* GetVectorUniquePointor() {
			return HKFormatError::ThrowErrorVersion<
			    std::unique_ptr<std::vector<std::unique_ptr<HKGeometryPMT> > >*>(__PRETTY_FUNCTION__, NULL);
		}

		//! Return the pointor to the vector
		virtual std::vector<std::unique_ptr<HKGeometryPMT> >* GetVectorPointor() {
			return HKFormatError::ThrowErrorVersion<std::vector<std::unique_ptr<HKGeometryPMT> >*>(
			    __PRETTY_FUNCTION__,
			    NULL);
		}

		template<class> friend class HKBranchManagerCollection;
#endif

		ClassDef(HKGeometryPMTCollection, 1);  //!< Class definition for ROOT
};

/************************************************
 *  HKGeometryPMTCollectionTemplate definition
 *
 *  Collection class wrapping a vector of HKGeometryPMT
 ***********************************************/
template<class T> class HKGeometryPMTCollectionT : public HKGeometryPMTCollection {

		static_assert(std::is_base_of<HKGeometryPMT, T>::value, "T must inherit from HKGeometryPMT");

	public:

		HKGeometryPMTCollectionT(bool initialize = true) {
			if(initialize) {
#ifdef HK_USE_ROOT7
				contents = std::make_shared<std::vector<std::unique_ptr<T> > >();
#else
				contents = std::make_unique<std::vector<std::unique_ptr<HKGeometryPMT> > >();
#endif
				this->Reset();
			}
		}  //!< Creator

		//! Constructor. Initialize the collection with specified number of entries
		/*!
		    \param entries integer. Number of entries
		*/
		HKGeometryPMTCollectionT(size_t entries) {
			this->Reset();
			this->resize(entries);
		}

		~HKGeometryPMTCollectionT() override { this->Reset(); }  //!< Destructor

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
		HKGeometryPMT* operator[](int i) override { return (*contents)[i].get(); }  // Use at(i) here too?

		//! Wrapper for std::vector at()
		/*!
		    \param i integer. Entry index
		*/
		HKGeometryPMT* at(int i) override { return contents->at(i).get(); }

		const HKGeometryPMT* at(int i) const { return contents->at(i).get(); }

		//! Wrapper for std::vector back()
		HKGeometryPMT* back() override { return contents->back().get(); }

		//! Wrapper for std::vector size()
		size_t size() const override { return contents->size(); }

		//! Wrapper for std::vector resize()
		/*!
		    \param entries integer. Number of entries
		*/
		void resize(size_t entries) override { 
			contents->resize(entries); 
#ifdef HK_USE_ROOT7
			for (auto& ptr : *contents) { ptr = std::make_unique<T>(); } 
#else
			for (auto& ptr : *contents) { ptr = std::make_unique<HKGeometryPMT>(); } 
#endif
		}

		//! Wrapper for std::vector push_back()
		/*!
		    \param in HKGeometryPMT*. Object to be inserted
		*/
		void push_back(HKGeometryPMT* in) override { contents->push_back(std::make_unique<T>(in)); }

		//! Wrapper for std::vector emplace_back()
		/*!
		    \param in HKGeometryPMT*. Object to be inserted
		*/
		void emplace_back(HKGeometryPMT* in) override { contents->emplace_back(std::make_unique<T>(in)); }

		//! Emplace back an empty HKGeometryPMT instance in the collection
		void AddGeometryPMT() override { contents->emplace_back(std::make_unique<T>()); }

	private:

#ifdef HK_USE_ROOT7
		std::shared_ptr<std::vector<std::unique_ptr<T> > > contents;  //!< Vector of HKGeometryPMT

		template<class, class, class> friend class HKDataMemberCollectionCreator;
#else
		//! Return the pointor to the unique_ptr
		std::unique_ptr<std::vector<std::unique_ptr<HKGeometryPMT> > >* GetVectorUniquePointor() override {
			return &contents;
		}

		//! Return the pointor to the vector
		std::vector<std::unique_ptr<HKGeometryPMT> >* GetVectorPointor() override { return contents.get(); }

		std::unique_ptr<std::vector<std::unique_ptr<HKGeometryPMT> > > contents;  //!< Vector of HKGeometryPMT

		template<class> friend class HKBranchManagerCollection;
#endif

		ClassDefOverride(HKGeometryPMTCollectionT, 1);  //!< Class definition for ROOT
};

#endif
