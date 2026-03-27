/*****************************************************************************************************/
//	LEAF.hpp																					
//	Authors: 	Guillaume Pronost (2019 - active) (pronost@km.icrr.u-tokyo.ac.jp)			
//				Benjamin Quilain (2019 - active)											
//				Lorenzo Perisse (2025 - active)												
//				Nicolas Moreau (???? - 2025)												
//				Shota Izumiyama (2020 - 2020)												
//	Original author: Benjamin Quilain														
//	Date: December 18th 2019																
//	Desc: Low-E Fitter for Hyper-K															
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

#ifndef LEAF_HPP
#define LEAF_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread> 
#include <mutex>

//* ROOT Headers
#include <TFile.h>
#include <TFitter.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TObject.h>
#include <TRandom3.h>
#include <TSpline.h>
#include <TPaletteAxis.h>
#include <TStopwatch.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>


#include "Constants/mPMTs.hpp"
#include "Enums/PMTType.hpp"
#include "HKGeometry.hpp"
#include "HKDarkNoise.hpp"
#include "HKHit.hpp"

#include "LEAFConfig.hpp"
#include "LEAFUtilities.hpp"
#include "SplineUtilities.hpp"

std::mutex mtx;

namespace LEAFStructures {

		struct FitterOutput {

			void Initialize() {
				vtx = ROOT::Math::XYZTVector();
				vtx_nll = 0.0;
				nll_r = 0.0;
				dir = ROOT::Math::XYZVector();
				dir_pre_fit = ROOT::Math::XYZVector();
				dir_nll = 0.0;
				energy = 0.0;
				total_charge = 0.0;
			}

			ROOT::Math::XYZTVector vtx; // leaf vertex
			double vtx_nll; // Vertex Likelihood
			double nll_r; // Goodness of fit

			ROOT::Math::XYZVector dir; // leaf Direction
			ROOT::Math::XYZVector dir_pre_fit; // leaf Direction pre-fit
			ROOT::Math::XYZVector dir_quick; // leaf Direction quick fit
			double dir_nll; // Direction Likelihood

			double energy; // leaf energy
			double total_charge; // raw total charge before affine correction
		};

		struct FitterPerformances {
			// Computations times
			double leaf_ct = 0.0; //Total leaf computational time
			double vtx_search_ct = 0.0; //Vertex coarse grid search computational time
			double vtx_minimize_ct = 0.0; //Vertex MINUIT minimisation computational time
			double dir_quick_search_ct = 0.0; //Direction using unit vectors from vector to PMTT - computational time
			double dir_search_ct = 0.0; //Direction coarse grid search computational time
			double dir_minimize_ct = 0.0; //Direction MINUIT minimisation computational time
			double energy_fit_ct = 0.0; //Energy finder computational time
		};	

		struct VtxCandidate {
			ROOT::Math::XYZTVector vtx;
			double nll = 0.0;
		};
}

namespace LEAFThreads {
	void SearchVertex_thread(int iStart, int iIte, int tolerance=1, bool likelihood=false, double lowerLimit=LEAFConfig::fSTimePDFLimitsQueueNegative, double upperLimit=LEAFConfig::fSTimePDFLimitsQueuePositive, int directionality=false);
	void MinimizeVertex_thread(int iStart, int iIte, std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, double lowerLimit=LEAFConfig::fSTimePDFLimitsQueueNegative, double upperLimit=LEAFConfig::fSTimePDFLimitsQueuePositive, int directionality = true, float FilterThreshold = 180);
};

namespace LEAFLikelihoods {
    void MinuitDirNLL(int& nDim, double* gout, double& NLL, double* par, int flg);
    void MinuitLikelihood(int& nDim, double * gout, double & NLL, double par[], int flg);
    void MinuitJointNLL(int& nDim, double * gout, double & NLL, double par[], int flg);
};

class LEAF 
{
	public:
		static LEAF*	GetME();
		void			DeleteME();

		/*****************/
  		/* CONFIGURATION */
		/*****************/

#ifdef HK_USE_ROOT7
		void Initialize(std::shared_ptr<HKGeometry const> lGeometry, std::shared_ptr<HKDarkNoise const> lDarkNoise, std::shared_ptr<HKGeometryPMTCollection const> lGeoPMT_ID, std::shared_ptr<HKGeometryPMTCollection const> lGeoPMT_mPMT);
#else
		void Initialize(const HKGeometry *lGeometry, const HKDarkNoise *lDarkNoise, const HKGeometryPMTCollection *lGeoPMT_ID, const HKGeometryPMTCollection *lGeoPMT_mPMT);
#endif
		void SetNThread(int iThread=N_THREAD);

		//* True Vertex isn't used for the fit, just to get feedback on how well the fit is at different steps
		void SetTrueVertexInfo(ROOT::Math::XYZTVector vtx);
		void SetTrueDirInfo(ROOT::Math::XYZVector dir);
		
		/****************/
  		/* MAIN PROCESS */
		/****************/

#ifdef HK_USE_ROOT7
		void LoadHitsCollection(std::shared_ptr<HKHitsCollection const> lHitCol_ID, std::shared_ptr<HKHitsCollection const> lHitCol_mPMT = nullptr);
#else
		void LoadHitsCollection(const HKHitsCollection* lHitCol_ID, const HKHitsCollection* lHitCol_mPMT = nullptr);
#endif

		//* Fitter Main Method. Process the whole fit. 
		LEAFStructures::FitterOutput MakeSequentialFit();

		//* slower and less performant for now
		LEAFStructures::FitterOutput MakeJointFit();

		//* Return performance information
		LEAFStructures::FitterPerformances GetPerformances() { return fPerformances; }

		/*****************/
  		/* VERTEX FITTER */
		/*****************/

		void SearchVertex_thread(int iStart, int iIte, int tolerance=1, bool likelihood=false, double lowerLimit=LEAFConfig::fSTimePDFLimitsQueueNegative, double upperLimit=LEAFConfig::fSTimePDFLimitsQueuePositive, int directionality=false);
		void MinimizeVertex_thread(int iStart, int iIte, std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, double lowerLimit=LEAFConfig::fSTimePDFLimitsQueueNegative, double upperLimit=LEAFConfig::fSTimePDFLimitsQueuePositive, int directionality = true, float FilterThreshold = 180);

		/***************/
  		/* LIKELIHOODS */
		/***************/

		//* Calculate likelihood without using PDF, but just using hits within a given timing window (hits "in-time").
		double Vertex_Score(ROOT::Math::XYZTVector vertexPosition, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);
		
		//* Contain the two functions above in one function, where usage of PDF or not can be set through the flag: likelihood = true/false
		double FindNLL(ROOT::Math::XYZTVector vertexPosition, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges=false, bool scaleDR=false, int directionality=false);

		double Vertex_Time_NLL(ROOT::Math::XYZTVector vertexPosition, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);

		double Dir_NLL(const ROOT::Math::XYZTVector& vertexPosition, const double theta_track, const double phi_track);
		
		double ComputeDirNLL_NoPDF(const ROOT::Math::XYZTVector& vertexPosition, const ROOT::Math::XYZVector& vertexDirection);

		double AngleNLL(ROOT::Math::XYZTVector vertexPosition, ROOT::Math::XYZVector vertexDirection);

		double FindNLLDirectionality(ROOT::Math::XYZTVector vertexPosition, int verbose, double lowerLimit, double upperLimit);

		double GoodnessOfFit(ROOT::Math::XYZTVector vertexPosition, double NLL_theory, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality);

	private:

		/*******************/
  		/* INTERNAL STRUCT */
		/*******************/

		struct EventInfo {
			int nHits = 0;
			double SignalOverNoise = 0.0;
			double NoiseIntegral = 0.0;
			double SignalIntegral = 0.0;
		};

		struct DirectionCandidate {
			double theta = 0.0;
			double phi = 0.0;
			double nll = 0.0;
		};

		struct JointFitCandidate {
			ROOT::Math::XYZTVector vtx;
			double theta = 0.0;
			double phi = 0.0;
			double nll = 0.0;
		};
		
		template <typename T>
		struct SortingNLL {
			bool operator()(T const &a, T const &b) const 
			{
				return a.nll < b.nll;
			}
		};

	private:
		
		LEAF();
		~LEAF();
		
		static LEAF* myFitter;
		std::vector<LEAFStructures::VtxCandidate> fVtxThreadOutput;
		std::vector<DirectionCandidate> fDirThreadOutput;

		/******************/
  		/* INITIALIZATION */
		/******************/

		void InitInputs();
		void MakePositionList();	
		void LoadSplines();

		/****************/
  		/* DATA PROCESS */
		/****************/

		void GenerateEventInfo(double lowerLimit, double upperLimit);
		double CalculateSNR(ROOT::Math::XYZTVector vertex, double lowerLimit, double upperLimit);

		/*****************/
  		/* VERTEX FITTER */
		/*****************/

		void FitVertex(float FilterThreshold = 180);
		
		//? Coarse GRID Search Parallelized version
		std::vector<LEAFStructures::VtxCandidate> SearchVertex_Main(int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality);
		
		//? Second Step Parallelized version
		std::vector<LEAFStructures::VtxCandidate> MinimizeVertex_Main(std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates = 1, int tolerance = 1, int verbose=0, bool likelihood=false, bool average=false, double lowerLimit=LEAFConfig::fSTimePDFLimitsQueueNegative, double upperLimit=LEAFConfig::fSTimePDFLimitsQueuePositive, int directionality = true, float FilterThreshold = 180);

		/********************/
  		/* DIRECTION FITTER */
		/********************/

		// Global fitter function that selects canidates and call the minimize method
		ROOT::Math::XYZVector FitDirection(const ROOT::Math::XYZTVector& vertex, bool searchPrior = false);

  		// Determine direction using the sum of all unit vectors from vertex to PMTs. Can be used as a prefit
		ROOT::Math::XYZVector FitDirectionQuick(const ROOT::Math::XYZTVector& vertex);

  		// Coarse GRID Search
		std::vector<DirectionCandidate> SearchDirection(const ROOT::Math::XYZTVector& vertex, int tolerance);
		
  		// LLH Minimization Search
		DirectionCandidate MinimizeDirection(const ROOT::Math::XYZTVector& vertex, std::vector<DirectionCandidate>& candidates, ROOT::Math::XYZTVector limits, int verbose);


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

	private:

		/*********************/
  		/* Private variables */
		/*********************/
				
		LEAFStructures::FitterOutput fOutput;
		LEAFStructures::FitterPerformances fPerformances; // additional infos about the fit
		ROOT::Math::XYZVector fDirectionFilter;
		bool fDirectionFilterSet;
		
		//! Number of threads to use for parallelization
		int fThread;

		//! Random number generator
		TRandom3 * fRand;

		//! Dark rate for each PMT type, in Hz
		std::array<double, static_cast<size_t>(PMTType::kNumPMTTypes)> fDarkRate_ns;

		//! List of position for grid coarse search
		std::vector<ROOT::Math::XYZTVector> fPositionList;

		//! True vertex position, use for debug
		ROOT::Math::XYZTVector fTrueVtxPos;

		//! True vertex direction, use for debug
		ROOT::Math::XYZVector fTrueDir;

		//! Splines
		std::array<TSpline3*, static_cast<size_t>(PMTType::kNumPMTTypes)> fSplineTimePDFConv;
		std::array<TSpline3*, static_cast<size_t>(PMTType::kNumPMTTypes)> fSplineTimePDFQueue;
		std::array<TSpline3*, static_cast<size_t>(PMTType::kNumPMTTypes)> fSplineTimePDFDarkRate;
		std::array<TSpline3*, static_cast<size_t>(PMTType::kNumPMTTypes)> fDirectionPDF;
		std::array<std::array<TGraph2D*,static_cast<size_t>(mPMTConstants::kNGroups)>, static_cast<size_t>(PMTType::kNumPMTTypes)> gPMTDirectionality_2D;
		std::array<TF1*, static_cast<size_t>(PMTType::kNumPMTTypes)> fDistResponsePMT;

		std::array<std::array<double, static_cast<size_t>(mPMTConstants::kNGroups)>, static_cast<size_t>(PMTType::kNumPMTTypes)> fDarkRate_dir_proba;
		double fLastLowerLimit;
		double fLastUpperLimit;

		//! 
		std::array<EventInfo, static_cast<size_t>(PMTType::kNumPMTTypes)> fEventInfo;

		//! DataModel
#ifdef HK_USE_ROOT7
		//! Geometry object
		std::shared_ptr<HKGeometry const> fGeometry;

		//! Dark noise object
		std::shared_ptr<HKDarkNoise const> fDarkNoise;

		//! PMT Informations
		std::array<std::shared_ptr<HKGeometryPMTCollection const>, static_cast<size_t>(PMTType::kNumPMTTypes)> fGeoPMTs;

		//! Hit Collections
		std::array<std::shared_ptr<HKHitsCollection const>, static_cast<size_t>(PMTType::kNumPMTTypes)> fHitsCollection;
		
#else
		//! Geometry object
		const HKGeometry* fGeometry;

		//! Dark noise object
		const HKDarkNoise* fDarkNoise;

		//! PMT Informations
		std::array<const HKGeometryPMTCollection*, static_cast<size_t>(PMTType::kNumPMTTypes)> fGeoPMTs;

		//! Hit Collections
		std::array<const HKHitsCollection*, static_cast<size_t>(PMTType::kNumPMTTypes)> fHitsCollection;
#endif
};

#endif
