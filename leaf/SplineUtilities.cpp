#include "SplineUtilities.hpp"

double SplineUtilities::SplineIntegral(TSpline3 *s, double start, double end, double stepSize) {
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * s->Eval(i + stepSize / 2);
	}
	return integral;
}

double SplineUtilities::SplineIntegralAndSubstract(TSpline3 *s0, TSpline3 *s1, double start, double end, double stepSize) {
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * (s0->Eval(i + stepSize / 2) - s1->Eval(i + stepSize / 2));
	}
	return integral;
}

double SplineUtilities::SplineIntegralExpo(TSpline3 *s, double start, double end, double sigma, double stepSize) {
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * s->Eval(i + stepSize / 2) * TMath::Gaus(i + stepSize / 2, 0, sigma);
	}
	return integral;
}