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
    int stepOneContainsTrueVtx;
    int stepOneContainsTrueDir;
    std::vector<double> Vtx;
    double NLL;
    double DNLL;
    double NLLR;

    int InTime;
    
    double True_NLLDiff;
    double True_TimeDiff;
    double True_TistDiff;

    double Energy;
    double TotalCharge;
    std::vector<double> Dir;
    std::vector<double> SNRList; 

    std::vector<double> MyDir;
    std::vector<double> Quick_Dir;

    std::vector<double> VtxHitHangles;
    std::vector<double> AngleCorrections;

    double Leaf_ComputeTime;
    double Vtx_Search_ComputeTime;
    double Vtx_Minimize_ComputeTime;
    double Dir_Quick_Search_ComputeTime;
    double Dir_Search_ComputeTime;
    double Dir_Minimize_ComputeTime;
    double Energy_Fit_ComputeTime;
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
    // std::vector<double> Vtx;
    // std::vector<double> Dir;
    double NLL;
};

struct FitPosition 
{
    std::vector<double> Vtx;
    double NLL;
};