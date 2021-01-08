#include "HitCollection.h"

template<typename T>
HitCollection<T>::HitCollection() {
	timestamp = 0;
	timestamp_last = 0;
	first_unique = 0;
}

template<typename T>
HitCollection<T>::~HitCollection() {
	this->Clean();
}

template<typename T>
void HitCollection<T>::Clean() {
	timestamp = 0;
	timestamp_last = 0;
	first_unique = 0;
	hits.clear();	
}

template<typename T>
void HitCollection<T>::SortByTime() {
	std::sort(hits.begin(), hits.end(), typename T::SortFunctor_HitTime());
}

template<typename T>
void HitCollection<T>::SortByTimeOfFlight() {
	std::sort(hits.begin(), hits.end(), typename T::SortFunctor_HitTimeOfFlight());
}

template<>
void HitCollection<Hit>::SortByTimeOfFlight() {
	//do nothing
	//std::cout << "FATAL: Call HitCollection<T>::SortByTimeOfFlight() for type: " << typeid(T).name() << std::endl;
}


template<typename T>
void HitCollection<T>::Append(const HitCollection<T> lHC) {
	// Adapted from SubSample append by T.Dealtry

	if ( hits.size() == 0 ) {
		// First time HitCollection is filled
		hits = lHC.hits;
		timestamp = lHC.timestamp;
		timestamp_last = lHC.timestamp_last;
		first_unique = lHC.first_unique;
	}
	else {
		// Need to shift the hit times by the difference of timestamp offsets
		TimeDelta::short_time_t time_shift = (lHC.timestamp - timestamp) / TimeDelta::ns;
	
		for ( unsigned int i = 0; i < lHC.Size(); i++ ) {
			T lHit = lHC.hits[i];
			lHit.T += time_shift;
		
			hits.push_back( lHit );
		}
		
		this->SortByTime();
		
		timestamp_last = timestamp + TimeDelta(hits.back().T);
	}
}

template<typename T>
void HitCollection<T>::ResetTimes() {

	// Check if there is any hit
	if ( hits.size() == 0 ) {
		//std::cout << "Empty collection" << std::endl;
		return;
	}
	
	this->SortByTime();
	
	// Get First Hit time
	TimeDelta::short_time_t time_shift = hits[0].T;
	
	// Check if time already reset
	if ( time_shift == 0 ) {
		//std::cout << "Sorted collection" << std::endl;
		// Reset only timestamp_last:
		timestamp_last = timestamp + TimeDelta(hits.back().T);
		return;
	}
	
	timestamp += time_shift;		
	
	for ( unsigned int i = 0; i < hits.size(); i++ ) {
		hits[i].T -= time_shift;		
	}
	
	timestamp_last = timestamp + TimeDelta(hits.back().T);
}

template<typename T> 
void HitCollection<T>::CopyCollection(const HitCollection<Hit> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.clear();
	
	hits.assign( lHC.hits.begin(), lHC.hits.end() );
}

template<typename T> 
void HitCollection<T>::CopyCollection(const HitCollection<HitExtended> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.clear();
	
	hits.assign( lHC.hits.begin(), lHC.hits.end() );
}

template<typename T> 
void HitCollection<T>::CopyDAQ(const HitCollection<T> lHC, const HKAA::DAQType lDAQ) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.clear();
	
	// copy:
	switch ( lDAQ ) {
		case HKAA::kDAQNone:
			std::cout << "[ERROR] CopyDAQ: No DAQ selected ! (" << lDAQ << ")" << std::endl;
			break;
		case HKAA::kIDDAQ_BnL:
			//std::cout << "[NOTICE] CopyDAQ: kIDDAQ_BnL (" << lDAQ << ")" << std::endl;
			std::copy_if ( lHC.hits.begin(), lHC.hits.end(), std::back_inserter(hits), [](T tHit) { return Astro_GetIfBnLPMT(tHit.PMT); } );
			break;
		case HKAA::kIDDAQ_3inch:
			//std::cout << "[NOTICE] CopyDAQ: kIDDAQ_3inch (" << lDAQ << ")" << std::endl;
			std::copy_if ( lHC.hits.begin(), lHC.hits.end(), std::back_inserter(hits), [](T tHit) { return Astro_GetIfmPMT(tHit.PMT); } );
			break;
		case HKAA::kIDDAQ_Hybrid:
			//std::cout << "[NOTICE] CopyDAQ: kIDDAQ_Hybrid (" << lDAQ << ")" << std::endl;
			std::copy_if ( lHC.hits.begin(), lHC.hits.end(), std::back_inserter(hits), [](T tHit) { return !Astro_GetIfODPMT(tHit.PMT); } );
			break;
		case HKAA::kODDAQ:
			//std::cout << "[NOTICE] CopyDAQ: kODDAQ (" << lDAQ << ")" << std::endl;
			std::copy_if ( lHC.hits.begin(), lHC.hits.end(), std::back_inserter(hits), [](T tHit) { return Astro_GetIfODPMT(tHit.PMT); } );
			break;
		default:
			std::cout << "[FATAL] CopyDAQ: this DAQ level is unknown ! (" << lDAQ << ")" << std::endl;
			exit(0);
	}	
}


template<typename T> 
void HitCollection<T>::FilterTime(const HitCollection<T> lHC, const TimeDelta tStart, const TimeDelta tEnd) {

	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.clear();
	
	for ( unsigned int i = 0; i < lHC.hits.size(); i++ ) {
	
		T tHit = lHC.hits[i];
		TimeDelta tTime = timestamp + TimeDelta(tHit.T);
		
		
		if ( tTime >= tStart && tTime <= tEnd ) 
			hits.push_back(tHit);		
	}
}

template<typename T> 
void HitCollection<T>::SetCollection(const HitCollection<T> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.clear();
	hits = lHC.hits;
}

template class HitCollection<Hit>;
template class HitCollection<HitExtended>;

