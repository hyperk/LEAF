/*****************************************************************************************************/
/**	LEAF.cc											**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)					**/
/**	Original author: Benjamin Quilain								**/
/**	Date: December 18th 2019									**/
/**	Desc: Low-E Fitter for Hyper-K								**/
/*****************************************************************************************************/

#include "LEAF.hpp"
#include <TStopwatch.h>

LEAF *LEAF::myFitter = NULL;

/*****************************************************************************************************/
/* GENERAL LEAF METHODS AND STRUCTURES FOR FITS */
/*****************************************************************************************************/

LEAF::LEAF() {
	myFitter = this;
	fDirectionFilterSet = false;
}

LEAF::~LEAF() {
	delete fRand;
	myFitter = NULL;
}

LEAF *LEAF::GetME() {
	if (myFitter)
		return myFitter;

	myFitter = new LEAF();
	return myFitter;
}

void LEAF::DeleteME() {
	if (myFitter) delete myFitter;
}

LEAFStructures::FitterOutput LEAF::MakeSequentialFit() {
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	fOutput.Initialize();
	
	if(fHitsCollection[int(PMTType::kID)]->size() == 0 && fHitsCollection[int(PMTType::kmPMT)] == 0) return fOutput;
	
	//* Vertex
	this->FitVertex();

	//* Goodness of fit
	fOutput.nll_r = this->GoodnessOfFit(fOutput.vtx, fOutput.vtx_nll, LEAFConfig::fMinimizeLimitsNegative, LEAFConfig::fMinimizeLimitsPositive, true, false, LEAFConfig::fUseDirectionality); //* Goodness of Fit
	
	//* Direction

	//? three direction methods are stored, one can use only the Fit Direction method
	// FitDirectionQuick(fOutput.Vtx); // prefit

	// // the prefit is stored in fOutput.Dir
	// FitDirection(fOutput.Vtx, fHitCollection->Size(), true); // optimization with one candidate that is the prefit
	// fOutput.MyDir = fOutput.Dir;
	// FitDirection(fOutput.Vtx, fHitCollection->Size(), false); // optimization without prefit, 30 candidates used

	this->FitDirectionQuick(fTrueVtxPos); // prefit

	// the prefit is stored in fOutput.Dir
	this->FitDirection(fTrueVtxPos, true); // optimization with one candidate that is the prefit
	fOutput.dir_pre_fit = fOutput.dir;
	this->FitDirection(fTrueVtxPos, false); // optimization without prefit, 30 candidates used
	

	//* Energy
	FitEnergy();

	timer.Stop();
	fPerformances.leaf_ct = timer.RealTime();

	return fOutput;
}

LEAFStructures::FitterOutput LEAF::MakeJointFit() {
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	fOutput.Initialize();

	std::vector<JointFitCandidate> Candidates = SearchVertexAndDir();
	JointFitCandidate finalCandidate = MinimizeVertexAndDir(Candidates);

	fOutput.vtx = finalCandidate.vtx;
	fOutput.dir = ROOT::Math::XYZVector(ROOT::Math::Polar3DVector(1.,finalCandidate.theta, finalCandidate.phi).Unit() );
	
	fOutput.vtx_nll = finalCandidate.nll;
	fOutput.dir_nll = this->Dir_NLL(fOutput.vtx, finalCandidate.theta, finalCandidate.phi);
	fOutput.nll_r = this->GoodnessOfFit(fOutput.vtx, fOutput.vtx_nll, LEAFConfig::fMinimizeLimitsNegative, LEAFConfig::fMinimizeLimitsPositive, true, false, LEAFConfig::fUseDirectionality);
	
	//* Energy, alone
	this->FitEnergy();

	timer.Stop();
	fPerformances.leaf_ct = timer.RealTime();

	return fOutput;
}

// Functions implementation was divided among the following files to make the code more readable
#include "LEAFSrcInputs.cpp"
#include "LEAFSrcFitter.cpp"
#include "LEAFSrcLikelihoods.cpp"