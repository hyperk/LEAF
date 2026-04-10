#ifndef LEAFUTILITIES_HPP
#define LEAFUTILITIES_HPP

#include <iostream>
#include <numeric>
#include <vector>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

#ifndef C_WATER_CM_NS
#define C_WATER_CM_NS   21.849964 // light celerity in water, in centimeter per ns.
#endif

namespace LEAFUtilities {
	
	double GetDistance(ROOT::Math::XYZVector a, ROOT::Math::XYZVector b);
	double GetLength(ROOT::Math::XYZVector v);
	
	double GetTimeOfFlight(ROOT::Math::XYZVector a, ROOT::Math::XYZVector b);
	double GetTimeOfFlight(ROOT::Math::XYZVector v);
	double GetTimeOfFlight(double distance);

	double GetResidual(ROOT::Math::XYZVector pmtPos, ROOT::Math::XYZVector vertexPos, double hitTime, double vertexTime);
	double GetResidual(ROOT::Math::XYZVector distToPMT, double hitTime, double vertexTime);
	double GetResidual(double tof, double hitTime, double vertexTime);

	double EstimateBandeWidth(const std::vector<double>& residuals); 

	double GaussianKernel(double x, double bandwidth);
	double KDE_Estimate(const std::vector<double>& residuals, double t_i, double bandwidth);

}; // namespace LEAFUtilities

#endif