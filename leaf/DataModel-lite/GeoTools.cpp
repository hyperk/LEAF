#include "GeoTools.h"

void GeoTools::CrossProduct( double* a, double* b, double* c ) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
} 

void GeoTools::Normalize( double* a ) {

	double dL = Astro_GetLength(a);
	a[0] /= dL;
	a[1] /= dL;
	a[2] /= dL;
} 
