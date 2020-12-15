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
			
			timestamp_last = timestamp + TimeDelta(lHit.T);
		}
	}
}

template<typename T> 
void HitCollection<T>::CopyCollection(const HitCollection<Hit> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.assign( lHC.hits.begin(), lHC.hits.end() );
}

template<typename T> 
void HitCollection<T>::CopyCollection(const HitCollection<HitExtended> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits.assign( lHC.hits.begin(), lHC.hits.end() );
}

template<typename T> 
void HitCollection<T>::SetCollection(const HitCollection<T> lHC) {
		
	timestamp = lHC.timestamp;
	timestamp_last = lHC.timestamp_last;
	first_unique = lHC.first_unique;
	
	hits = lHC.hits;
}

template class HitCollection<Hit>;
template class HitCollection<HitExtended>;

