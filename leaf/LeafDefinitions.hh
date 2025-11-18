#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread> 
#include <mutex>

//* WCSim Headers
#include "WCSimRootGeom.hh"

//* ROOT Headers
#include "TFile.h"
#include "TFitter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TPaletteAxis.h"

//* DataModel informations
#include "Geometry.h"
#include "HitCollection.h"

struct FitterOutput 
{
    std::vector<double> Vtx; // leaf vertex
    double NLL; // Vertex Likelihood
    double NLLR; // Goodness of fit

    std::vector<double> Dir; // leaf Direction
    double DNLL; // Direction Likelihood

    double Energy; // leaf energy
    double TotalCharge; // raw total charge before affine correction
    
    std::vector<double> MyDir;
    std::vector<double> Quick_Dir;
};

struct FitterOutputProps
{
    // Computations times
    double Leaf_ComputeTime;//Total leaf computational time
    double Vtx_Search_ComputeTime;//Vertex coarse grid search computational time
    double Vtx_Minimize_ComputeTime;//Vertex MINUIT minimisation computational time
    double Dir_Quick_Search_ComputeTime;//Direction using unit vectors from vector to PMTT - computational time
    double Dir_Search_ComputeTime;//Direction coarse grid search computational time
    double Dir_Minimize_ComputeTime;//Direction MINUIT minimisation computational time
    double Energy_Fit_ComputeTime;//Energy finder computational time
};

struct DirectionCandidate 
{
    double theta;
    double phi;
    double DNLL;
};

struct VtxCandidate 
{
    double X;
    double Y;
    double Z;
    double T;
    double NLL;
    double SNR;  // SNR value
};

struct JointFitCandidate
{
    VtxCandidate VtxPart;
    DirectionCandidate DirPart;
    double NLL;
};

struct FitPosition 
{
    std::vector<double> Vtx;
    double NLL;
};
