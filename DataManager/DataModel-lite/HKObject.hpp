#ifndef HK_OBJECT_HPP
#define HK_OBJECT_HPP

/************************************************
 *  HKObject.hpp
 *   - Author: Guillaume Pronost
 *   - Date: 2024/07/10
 ***********************************************/

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TObject.h>

#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#define UNUSED_PARAM __attribute__((unused))

namespace HKFormatError {

	//! Generic function to throw a runtime error when a virtual function is called.
	/*!
	    \param function_name string. Name of the virtual function
	    \param default_value type T. Default value, in case error is catched.
	    \return type T default_value defined in parameters
	*/
	template<typename T> T ThrowErrorVersion(std::string function_name, T default_value) {
		throw std::runtime_error(function_name + " doesn't exist in this format version");
		return default_value;
	}
}  // namespace HKFormatError

namespace HKFormatTools {

	//! Generic function to generate branch name
	/*!
	    \param suffix string. Suffix in the branch name
	    \return Branch name
	*/
	template<class T> std::string GetBranchName(std::string suffix = "") {

		std::stringstream stream;
		if(suffix == "") {
			stream << T::GetBaseNameStatic() << "V" << T::GetVersionStatic();
		}
		else {
			stream << T::GetBaseNameStatic() << "_" << suffix << "V" << T::GetVersionStatic();
		}
		return stream.str();
	}

	//! Generic function to generate registration name
	/*!
	    \param suffix string. Suffix in the branch name
	    \return Branch name
	*/
	template<class T> std::string GetRegistrationName(std::string suffix = "") {

		std::stringstream stream;
		if(suffix == "") {
			stream << T::GetBaseNameStatic();
		}
		else {
			stream << T::GetBaseNameStatic() << "_" << suffix;
		}
		return stream.str();
	}
}  // namespace HKFormatTools

/**********************************************************/

/************************************************
 *  HKObject class definition
 *
 *  Default class from which HK data format class derive.
 *
 *  Define the default functions each class should have
 ***********************************************/
class HKObject {
	public:

		HKObject() {}  //!< Creator

		virtual ~HKObject() = 0;  //!< Destructor

		// Each object needs a reset function
		virtual void Reset() = 0;  //!< Reset function, called after each event

		// Each object needs a version
		virtual unsigned int GetVersion() const = 0;  //!< Return the class version number

		// Each object needs a name
		virtual std::string GetBaseName() const = 0;  //!< Return the class base name

		std::string GetBranchName(std::string suffix = "") const {
			std::stringstream stream;
			if(suffix == "") {
				stream << this->GetBaseName() << "V" << this->GetVersion();
			}
			else {
				stream << this->GetBaseName() << "_" << suffix << "V" << this->GetVersion();
			}
			return stream.str();
		}  //!< Return the branch name to be used in input/output

		static constexpr bool isCollection() {
			return false;
		}  //!< Static function. Return if the class is a collection (i.e. a vector of object)

		ClassDef(HKObject, 1);  //!< Class definition for ROOT
};

/************************************************
 *  HKObjectCollection class definition
 *
 *  Collection class wrapping a vector of HKObject
 *
 *  Default class from which HK data format collection class derive.
 *
 *  Define the default functions each class should have
 ***********************************************/
class HKObjectCollection : public HKObject {

	public:

		HKObjectCollection() {}  //!< Creator

		virtual ~HKObjectCollection()           = 0;  //!< Destructor
		virtual void Reset()                    = 0;  //!< Reset function, called after each event
		virtual unsigned int GetVersion() const = 0;  //!< Return the contained class version number

		static constexpr bool isCollection() {
			return true;
		}  //!< Static function. Return if the class is a collection (i.e. a vector of object)

		ClassDef(HKObjectCollection, 1);  //!< Class definition for ROOT
};

#endif
