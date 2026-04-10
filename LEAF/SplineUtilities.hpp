#ifndef SPLINE_UTILITIES_HPP
#define SPLINE_UTILITIES_HPP

#include <TMath.h>
#include <TSpline.h>

namespace SplineUtilities {
    double SplineIntegral(TSpline3 *s, double start, double end, double stepSize=5e-1);
    double SplineIntegralAndSubstract(TSpline3 *s0, TSpline3 *s1, double start, double end, double stepSize=5e-1);
    double SplineIntegralExpo(TSpline3 *s, double start, double end, double sigma, double stepSize=5e-1);
}; // namespace SplineUtilities

#endif