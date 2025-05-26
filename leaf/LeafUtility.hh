#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

//* DataModel informations
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafConfig.hh"
#include "LeafInputs.hh"
#include "LeafDefinitions.hh"

inline bool SortOutputVector ( const std::vector<double>& v1, const std::vector<double>& v2 ) 
{ 
	return v1[4] < v2[4]; 
}

double Distance3D(std::vector<double> point1, std::vector<double> point2);
		
void Normalize(double a[3]);

void Normalize(std::vector<double>& vector);

double dot(const std::vector<double>& A, const std::vector<double>& B);

double calculateDistance(const std::vector<double>& A, const std::vector<double>& B);

std::vector<double> PolarToCartesianNorm(std::vector<double> polarVector);

std::vector<double> CartesianToPolarNorm(const std::vector<double>& cartVector);

bool CanidatesMatch(std::vector<std::vector<double>>* candidates, std::vector<double> point, std::vector<double> fTrueVtxPos);

bool ContainsTrueVtx(std::vector<VtxCandidate>* candidates, std::vector<double> fTrueVtxPos);

bool ContainsTrueDir(std::vector<DirectionCandidate>* candidates, std::vector<double> fTrueDir);

//* Estimate P_data(t_i) using KDE
double KDE_Estimate(const std::vector<double>& residuals, double t_i, double bandwidth);

double EstimateBandeWidth(const std::vector<double>& residuals);

double ComputeResidualTime(std::vector<double> vertexPos, double originTime, Hit lHitt);

//* Gaussian kernel function
double GaussianKernel(double x, double bandwidth);

//* Signal noise ratio
VtxCandidate ComputeCandidateSNR(const std::vector<double>& vertex, double lowerLimit, double upperLimit);

void VectorVertexPMT(std::vector<double> vertex, int iPMT, double* dAngles );