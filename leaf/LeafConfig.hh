#pragma once

#define NPMT_CONFIGURATION 	1

// Hit definition:
#define NormalPMT		0 // B&L hit
#define MiniPMT		1 // 3" PMT hit
#define AllPMT			2

#define VERBOSE 		0

#undef VERBOSE_VTX // In SearchVertex
#undef VERBOSE_NLL // In Likelihood Computation

#define VTX_X			0
#define VTX_Y			1
#define VTX_Z			2
#define VTX_T			3

#define N_THREAD		12 // Default for sukap

#include "WCSimRootGeom.hh"
#include "Geometry.h"

// extern constexpr float fCVacuum = 3e8 * 1e2 / 1e9; // speed of light, in centimeter per ns.
// extern constexpr float fNIndex = 1.373;		
// extern constexpr double fLightSpeed = fCVacuum / fNIndex;

//* Geometry & Constants
extern double fTankRadius;
extern double fTankHeight;
extern double fTankHalfHeight;
extern double fLightSpeed;

//* PDFs
extern double fPDFNorm_fullTimeWindow;
extern double fSTimePDFLimitsQueueNegative;
extern double fSTimePDFLimitsQueuePositive;
extern double fSTimePDFLimitsQueueNegative_fullTimeWindow;
extern double fSTimePDFLimitsQueuePositive_fullTimeWindow;

//* Overall Fit Parameters
extern double fTimeWindowSizeFull;
extern double fIntegrationTimeWindow;
extern bool   fStepByStep;	
extern bool   fUseDirectionality;
extern bool	DoubleNLL;
extern bool   fLimit_mPMT;
// extern double KELimitsPos;
// extern double KELimitsNeg;
extern int    fAveraging;
extern bool   fHighEnergy;
extern double fSearchVtxStep;
extern double fSearchVtxTolerance;

//* Step 1 Specific Parameters
extern double fHitTimeLimitsNegative;
extern double fHitTimeLimitsPositive;

//* Step 2 Specific Parameters
extern double fMinimizeLimitsNegative;
extern double fMinimizeLimitsPositive;

//* Direction Fit Parameters
extern bool DirTakeAll;

void InitConfig(const Geometry *lGeometry);