#ifndef GeoTools_hh
#define GeoTools_hh

#include "Environments.h"

/**
* \namespace GeoTools
 *
 * This namespace hold some useful function for computation
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/18 16:26:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 
 
namespace GeoTools {

	void CrossProduct( double* a, double* b, double* c );
	void Normalize( double* a );
}

#endif
