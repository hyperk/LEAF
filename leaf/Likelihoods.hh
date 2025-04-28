#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>
#include <cmath>
#include <numeric>

#include "TGraph2D.h"

#include "WCSimRootGeom.hh"
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafInputs.hh"
#include "LeafSplines.hh"
#include "LeafUtility.hh"
#include "LeafDefinitions.hh"

class Likelihoods
{
	public:
    struct SortingNLL
    {
        bool operator()(FitPosition const &a, FitPosition const &b) const 
        {
            return a.NLL < b.NLL;
        }
    };

    //* Calculate likelihood without using PDF, but just using hits within a given timing window (hits "in-time").
    static double Vertex_Score(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);
    
    static double FindNLL_NoLikelihood_Energy(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);
    
    //* Contain the two functions above in one function, where usage of PDF or not can be set through the flag: likelihood = true/false
    static double FindNLL(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition,int nhits, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges=false, bool scaleDR=false, int directionality=false);

    static double Vertex_Time_NLL(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);

    static double Dir_NLL(const HitCollection<Hit>* lHitCol, const std::vector<double>& vertexPosition, const double theta_track, const double phi_track, int nhits);
	
    static double ComputeDirNLL_NoPDF(const HitCollection<Hit>* lHitCol, const std::vector<double>& vertexPosition, const std::vector<double>& vertexDirection, int nhits);

    static double FindNLLDirectionality(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, int verbose, double lowerLimit, double upperLimit);

    static double GoodnessOfFit(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);

    //* Estimate P_data(t_i) using KDE
    static double KDE_Estimate(const std::vector<double>& residuals, double t_i, double bandwidth);

    static double EstimateBandeWidth(const std::vector<double>& residuals);
    
    private:

    static double ComputeResidualTime(std::vector<double> vertexPos, double originTime, Hit lHitt);

    //* Gaussian kernel function
    static double GaussianKernel(double x, double bandwidth);

};