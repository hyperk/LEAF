#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>

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

//WCSim Headers
#include "WCSimRootGeom.hh"
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafInputs.hh"


struct EventInfo 
{
    int hits;
    double SignaloverNoise;
    double NoiseIntegral;
    double SignalIntegral;
};

// extern bool fHighEnergySplines;

extern TSpline3 *	fSplineTimePDFQueue[NPMT_CONFIGURATION];
extern TSpline3 *	fSplineTimePDFDarkRate[NPMT_CONFIGURATION];
extern TSpline3* fDirectionPDF;
extern TSpline3* fHitAnglePDF;

// Histo
extern TGraph2D * 	gPMTDirectionality_2D[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];

// TF1
extern TF1 * 		fDistResponsePMT[NPMT_CONFIGURATION];

extern double fDarkRate_dir_proba[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];

extern double fLastLowerLimit;
extern double fLastUpperLimit;
// double fTimeWindowSizeFull;

void LoadSplines();

double SplineIntegral			(TSpline3 * s, 			double start, double end, 		double stepSize=5e-1);
double SplineIntegralAndSubstract	(TSpline3 * s0, TSpline3 * s1, 	double start, double end, 		double stepSize=5e-1);
double SplineIntegralExpo		(TSpline3 * s, 			double start, double end, double sigma, double stepSize=5e-1);
    
EventInfo MakeEventInfo(double lowerLimit, double upperLimit, int pmtType);