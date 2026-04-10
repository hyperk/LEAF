#ifndef LEAFCONFIG_HPP
#define LEAFCONFIG_HPP

#include <iostream>
#include <map>
#include <vector>
#include <mutex>
#include <TMath.h>
#include "Enums/PMTType.hpp"

#define VERBOSE 		0

#undef VERBOSE_NLL // In Likelihood Computation

#define N_THREAD		2 // Default for sukap
#ifndef C_WATER_CM_NS
#define C_WATER_CM_NS   21.849964 // light celerity in water, in centimeter per ns.
#endif


class LEAFConfig {

	public:
		static std::mutex mtx;
		
		//* Geometry & Constants
		static double fTankRadius;
		static double fTankHeight;
		static double fTankHalfHeight;
		static double fLightSpeed;

		//* PDFs
		static double fSTimePDFLimitsQueueNegative;
		static double fSTimePDFLimitsQueuePositive;
		static double fSTimePDFLimitsQueueNegative_fullTimeWindow;
		static double fSTimePDFLimitsQueuePositive_fullTimeWindow;
		static double fPDFNorm_fullTimeWindow;

		static double fDirectionPDF_maxResidual;
		static double fDirectionPDF_minResidual;

		//* Overall Fit Parameters
		static double fTimeWindowSizeFull;
		static double fIntegrationTimeWindow;
		static bool   fStepByStep;
		static bool   fUseDirectionality;
		static bool   fDoubleNLL;
		static bool   fLimit_mPMT;
		static int	fAveraging;
		static bool   fHighEnergy;
		static double fSearchVtxStep;
		static double fSearchVtxTolerance;

		//* Step 1 Specific Parameters
		static double fHitTimeLimitsNegative;
		static double fHitTimeLimitsPositive;

		//* Step 2 Specific Parameters
		static double fMinimizeLimitsNegative;
		static double fMinimizeLimitsPositive;

		//* Direction Fit Parameters
		static bool fDirTakeAll;
		static double fThetaStep; 
		static double fPhiStep; 
		static int fDirTolerance;

		//! List of active PMT types
		static std::vector<PMTType> fActivePMTTypes; 

		static void Initialize(double lTankRadius, double lTankHeight);

   
}; // class LEAFConfig

#endif