#include "LEAFUtilities.hpp"

double LEAFUtilities::GetDistance(ROOT::Math::XYZVector a, ROOT::Math::XYZVector b) {
	ROOT::Math::XYZVector distVector = (a - b);
	return sqrt(distVector.Dot(distVector));	
}

double LEAFUtilities::GetLength(ROOT::Math::XYZVector v) {
	return sqrt(v.Dot(v));	
}

double LEAFUtilities::GetTimeOfFlight(ROOT::Math::XYZVector a, ROOT::Math::XYZVector b) {
	return LEAFUtilities::GetDistance(a,b) / C_WATER_CM_NS;
}

double LEAFUtilities::GetTimeOfFlight(ROOT::Math::XYZVector v) {
	return LEAFUtilities::GetLength(v) / C_WATER_CM_NS;
}

double LEAFUtilities::GetTimeOfFlight(double distance) {
	return distance / C_WATER_CM_NS;
}

double LEAFUtilities::GetResidual(ROOT::Math::XYZVector pmtPos, ROOT::Math::XYZVector vertexPos, double hitTime, double vertexTime) {
	return hitTime - LEAFUtilities::GetTimeOfFlight(vertexPos, pmtPos) - vertexTime;
}

double LEAFUtilities::GetResidual(ROOT::Math::XYZVector distToPMT, double hitTime, double vertexTime) {
	return hitTime - LEAFUtilities::GetTimeOfFlight(distToPMT) - vertexTime;
}

double LEAFUtilities::GetResidual(double tof, double hitTime, double vertexTime) {
	return hitTime - tof - vertexTime;
}

double LEAFUtilities::EstimateBandeWidth(const std::vector<double>& residuals) {
	double sum = std::accumulate(residuals.begin(), residuals.end(), 0.0);
	double mean = sum / residuals.size();

	double sq_sum = std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0);
	double variance = sq_sum / residuals.size() - mean * mean;
	double stddev = std::sqrt(variance);
	return 1.06 * stddev * std::pow(residuals.size(), -1.0 / 5.0);
}

double LEAFUtilities::GaussianKernel(double x, double bandwidth) {
    return std::exp(-0.5 * x * x / (bandwidth * bandwidth)) / (bandwidth * std::sqrt(2.0 * M_PI));
}

double LEAFUtilities::KDE_Estimate(const std::vector<double>& residuals, double t_i, double bandwidth) {
    double sum = 0.0;
    for (const double& t_j : residuals) sum += GaussianKernel(t_i - t_j, bandwidth);
    return sum / residuals.size();
}