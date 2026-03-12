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

//* LEAF Headers
#include "Likelihoods.hpp"
#include "LeafInputs.hpp"
#include "LeafConfig.hpp"
#include "LeafSplines.hpp"
#include "LeafUtility.hpp"
#include "LeafDefinitions.hpp"



class LEAF 
{
	public:
		static LEAF*		GetME();
		void			DeleteME();

		static FitterOutput fOutput;
		static FitterOutputProps fOutputProps; // additional infos about the fit
		
		void Initialize(const Geometry* lGeometry);
		
		//* Fitter Main Method. Process the whole fit. 
		struct FitterOutput MakeSequentialFit(const HitCollection<Hit>* lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT=true);

		//* slower and less performant for now
		struct FitterOutput MakeJointFit(const HitCollection<Hit>* lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT=true);


		/*****************/
  		/* VERTEX FITTER */
		/*****************/

		void SearchVertex_thread(int iStart, int iIte, int nhits, int tolerance=1, bool likelihood=false, double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality=false);
		void MinimizeVertex_thread(int iStart, int iIte, std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits, int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality = true, std::vector<double>* fDirection_Filter = nullptr, float FilterThreshold = 180);


	private:
		
		LEAF();
		~LEAF();
		struct FitterOutput NewOutuput();
		void LoadHitCollection(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT);
		static LEAF* myFitter;
		std::vector< std::vector<double> > fThreadOutput;

		/*****************/
  		/* VERTEX FITTER */
		/*****************/

		void FitVertex(std::vector<double>* fDirection_Filter = nullptr, float FilterThreshold = 180);
		
		//? Coarse GRID Search Parallelized version
		std::vector<VtxCandidate> SearchVertex_Main(int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality);
		
		//? Second Step Parallelized version
		std::vector< std::vector<double> > MinimizeVertex_Main(std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits, int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, double lowerLimit=fSTimePDFLimitsQueueNegative, double upperLimit=fSTimePDFLimitsQueuePositive, int directionality = true, std::vector<double>* fDirection_Filter = nullptr, float FilterThreshold = 180);


		/********************/
  		/* DIRECTION FITTER */
		/********************/

		// Global fitter function that selects canidates and call the minimize method
		std::vector<double> FitDirection(const std::vector<double>& fixedVertexPosition, int nhits, bool searchPrior = false);

  		// Determine direction using the sum of all unit vectors from vertex to PMTs. Can be used as a prefit
		std::vector<double> FitDirectionQuick(const std::vector<double>& fixedVertexPosition);

  		// Coarse GRID Search
		std::vector<DirectionCandidate> SearchDirection(const std::vector<double>& fixedVertexPosition, int tolerance);
		
  		// LLH Minimization Search
		DirectionCandidate MinimizeDirection(const std::vector<double>& fixedVertexPosition, std::vector<DirectionCandidate>& candidates, int nhits, double *limits, int verbose);


		/*****************/
  		/* ENERGY FITTER */
		/*****************/

  		void FitEnergy();


		/****************/
  		/* JOINT FITTER */
		/****************/

		// Search combined candidates
  		std::vector<JointFitCandidate> SearchVertexAndDir();

		// Minimize these candidates
		JointFitCandidate MinimizeVertexAndDir(std::vector<JointFitCandidate> initialCandidates);
};

#endif
