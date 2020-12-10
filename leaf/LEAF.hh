/*****************************************************************************************************/
/**	LEAF.hh											**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)					**/
/**	Original author: Benjamin Quilain								**/
/**	Date: December 18th 2019									**/
/**	Desc: Low-E Fitter for Hyper-K								**/
/**     + Add Low-E functions from Super-K								**/
/*****************************************************************************************************/

/*****************************************************************************************************/
// This class contains a vertex finder for Low-E WCSim-based detectors through different methods
// The vertex finder relies on a likelihood of timing residual for each PMT hit.
// Time residuals are defined as:
// time - time-of-flight-assuming-straight-line-from-vertex-to PMT - time-of-vertex
/*****************************************************************************************************/
/*****************************************************************************************************/
// The standard method relies on two parts which are runned successively:
//
// 1. SearchVertex: A coarse grid search in the tank, where steps in space/time can be set by users.
// In default mode, this class does not use the full likelihood but just number of hits in time to 
// be faster
//
// 2. MinimizeVertex: A MINUIT-based search within a 4D sphere around a specific position and time
// Note that in default mode, the sphere radius should be set to be the same or similar to the step 
// size in SearchVertex. In default mode, the best candidates from SearchVertex are fed to 
// MinimizeVertex. The number of best candidates can be selected. // Both MinimizeVertex and 
// SearchVertex provides an output containing the best fits vertices, ordered from lower to higher 
// NLL. 
//
// Most of default parameters of method 1 and 2 can be set in LEAF.cc Init() method
/*****************************************************************************************************/


#ifndef LEAF_hh
#define LEAF_hh

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread> 
#include <mutex>

//WCSim Headers
#include "WCSimRootGeom.hh"

//ROOT Headers
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

//DataModel informations
#include "Geometry.h"
#include "HitCollection.h"

// Number of PMT configuration:
#define NPMT_CONFIGURATION 	2

// Hit definition:
#define NormalPMT		0 // B&L hit
#define MiniPMT		1 // 3" PMT hit
#define AllPMT			2

#define VERBOSE 		0

// Verbose level in functions:
#undef VERBOSE_VTX // In SearchVertex
#undef VERBOSE_NLL

#define VTX_X			0
#define VTX_Y			1
#define VTX_Z			2
#define VTX_T			3

#define N_THREAD		12 // Default for sukap

// Likelihood:
void MinuitLikelihood(int& nDim, double * gout, double & NLL, double par[], int flg);
void MinimizeVertex_CallThread(	int iStart, int iIte,	
					std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
					int nCandidates, int tolerance, int verbose, bool likelihood, bool average,
					double lowerLimit, double upperLimit, int directionality);

void SearchVertex_CallThread(	int iStart, int iIte,
				int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality);

inline bool SortOutputVector ( const std::vector<double>& v1, const std::vector<double>& v2 ) { 
	return v1[4] < v2[4]; 
} 
		
class LEAF {

	public:
		static LEAF*		GetME();
		void			DeleteME();
		
		void Initialize(const Geometry* lGeometry);
		
		void SetNThread(int iThread=N_THREAD) { fThread=iThread; }
		
		struct FitterOutput {
		
			double Vtx[4];
			double NLL;
			int InTime;
			
			double True_NLLDiff;
			double True_TimeDiff;
			double True_TistDiff;
		};
		

                // Fitter Main Method. Process the whole fit. 
		struct FitterOutput MakeFit(const HitCollection<Hit>* lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT=true);
		
		// NLL
                // Calculate likelihood function based on input PDF. For now, the function is based on time residuals.
		double FindNLL_Likelihood(std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);
                // Calculate likelihood without using PDF, but just using hits within a given timing window (hits "in-time").
		double FindNLL_NoLikelihood(std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);
                // Contain the two functions above in one function, where usage of PDF or not can be set through the flag: likelihood = true/false
		double FindNLL(		std::vector<double> vertexPosition,int nhits, bool likelihood, int verbose, double lowerLimit, double upperLimit, 
						bool killEdges=false, bool scaleDR=false, int directionality=false);
		
		void SearchVertex_thread(	int iStart, int iIte,	
						int nhits, 
						int tolerance=1, bool likelihood=false, 
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality=false);

		
		void MinimizeVertex_thread(	int iStart, int iIte,	
						std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
						int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false,
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality = true);
		
	private:
	
		LEAF();
		~LEAF();
		
		void Init();
		void LoadSplines();
		void LoadPMTInfo();
		void MakeMPMTReferencial(int iPMT);
		void MakePositionList();
		
		double GetDistanceOld(std::vector<double> A, std::vector<double> B);
		
		// Spline functions:
		double SplineIntegral			(TSpline3 * s, 			double start, double end, 		double stepSize=5e-1);
		double SplineIntegralAndSubstract	(TSpline3 * s0, TSpline3 * s1, 	double start, double end, 		double stepSize=5e-1);
		double SplineIntegralExpo		(TSpline3 * s, 			double start, double end, double sigma, double stepSize=5e-1);
			
		// NLL
		void MakeEventInfo(double lowerLimit, double upperLimit);

		double FindNLLDirectionality(std::vector<double> vertexPosition, int nhits, int verbose, double lowerLimit, double upperLimit);
		
		
		double FindDirectionTheta(std::vector<double> vertex,int tubeNumber, int verbose);

		// SearchVertex functions:
		std::vector< std::vector<double> > SearchVertex(		
						int nhits, 
						int tolerance=1, bool likelihood=false, 
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality=false);
						
		std::vector< std::vector<double> > SearchVertex_Main(		
						int nhits, 
						int tolerance=1, bool likelihood=false, 
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality=false);
						
						
		std::vector< std::vector<double> > SearchVertexFine(	
						std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,  
						int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, 
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality=false);
		
		// Minimize function:
		std::vector< std::vector<double> > MinimizeVertex(	
						std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
						int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false,
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality = true);
						
		std::vector< std::vector<double> > MinimizeVertex_Main(	
						std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
						int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false,
						double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality = true);

		
		void VectorVertexPMT( std::vector<double> vertex, int iPMT, double* dAngles );
	
	// Variables:
	private:
	
		struct FitPosition {
			std::vector<double> Vtx;
			double NLL;
		};
		struct SortingNLL { 
			bool operator()(FitPosition const &a, FitPosition const &b) const { 
				return a.NLL < b.NLL;
			}
		};
		
		struct EventInfo {
			int hits;
			double SignaloverNoise;
			double NoiseIntegral;
			double SignalIntegral;
		};
		
		static LEAF* myFitter;

		// Spline
		TSpline3 *	fSplineTimePDFQueue[NPMT_CONFIGURATION];
		TSpline3 *	fSplineTimePDFDarkRate[NPMT_CONFIGURATION];
		
		// Histo
		TGraph2D * 	gPMTDirectionality_2D[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];
				
		// TF1
		TF1 * 		fDistResponsePMT[NPMT_CONFIGURATION];
		
		
		double fDarkRate_dir_proba[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];
		
		static double fSTimePDFLimitsQueueNegative; 
		static double fSTimePDFLimitsQueuePositive; 
		double fSTimePDFLimitsQueueNegative_fullTimeWindow;
		double fSTimePDFLimitsQueuePositive_fullTimeWindow;
		double fTimeWindowSizeFull;
		
		double fMinimizeLimitsNegative;
		double fMinimizeLimitsPositive;
		
		double fHitTimeLimitsNegative;
		double fHitTimeLimitsPositive;
		
		
		// Inputs:
		const Geometry* fGeometry;
		const HitCollection<Hit>* fHitCollection;
		TimeDelta fTriggerTime;
		TimeDelta fTimeCorrection;
		
		// PMT Informations
		const std::vector<PMTInfo> *fPMTList;
		double fDarkRate_ns[NPMT_CONFIGURATION];
		
		// Detector geometry:
		double fTankRadius;
		double fTankHeight;
		double fTankHalfHeight;
		double fLightSpeed;
		
		// Fitter parameters:
		double fIntegrationTimeWindow;
		bool   fStepByStep;
		bool   fUseDirectionality;
		bool   fLimit_mPMT;
		int    fAveraging;
		bool   fHighEnergy;
		double fSearchVtxStep;
		double fSearchVtxTolerance;
		
		// Input/Output
		//std::vector<double> fTrueVtxPos;
		//std::vector< std::vector<double> > fTrueVtxPosDouble;
		double fPDFNorm_fullTimeWindow;
		std::vector< std::vector<double> > fRecoVtxPosFinal;
		
		
		
		// Random generator
		TRandom3 * fRand;					
		
		double fLastLowerLimit;
		double fLastUpperLimit;
		struct EventInfo fEventInfo[NPMT_CONFIGURATION];
		
		int fThread;
		
		std::vector< std::vector<double> > fPositionList;
		
		
		std::vector< std::vector<double> > fThreadOutput;
		
};

#endif
