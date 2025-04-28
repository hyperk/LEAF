#include "LeafConfig.hh"

//* Geometry & Constants
double fTankRadius;
double fTankHeight;
double fTankHalfHeight;
double fLightSpeed;

//* PDFs
double fSTimePDFLimitsQueueNegative = -3.;
double fSTimePDFLimitsQueuePositive = 4.;
double fSTimePDFLimitsQueueNegative_fullTimeWindow = 0.;
double fSTimePDFLimitsQueuePositive_fullTimeWindow = 0.;
double fPDFNorm_fullTimeWindow = 0.;

//* Overall Fit Parameters
double fTimeWindowSizeFull = 1500.;
double fIntegrationTimeWindow = 50.;
bool   fStepByStep = false;
bool   fUseDirectionality = false;
bool	DoubleNLL = false;
bool   fLimit_mPMT = true;
// double KELimitsPos;
// double KELimitsNeg;
int    fAveraging = 20;
bool   fHighEnergy = false;
double fSearchVtxStep = 300.;
double fSearchVtxTolerance = 60;

//* Step 1 Specific Parameters
double fHitTimeLimitsNegative = -5;
double fHitTimeLimitsPositive = -7;

//* Step 2 Specific Parameters
double fMinimizeLimitsNegative = -700;
double fMinimizeLimitsPositive = 1000;

//* Direction Fit Parameters
bool DirTakeAll = false;

void InitConfig(const Geometry *lGeometry)
{
    fTankRadius = lGeometry->detector_radius;
	fTankHeight = lGeometry->detector_length;
	fTankHalfHeight = fTankHeight / 2.;


    fStepByStep = false; // Step by step mode was not test and is not supported by multithreading
	fUseDirectionality = false;
	fHighEnergy = false;
	DoubleNLL = true;

	// If true, we apply same cuts for timing as B&L PMT when searching vertex.
	// Otherwise, we do not apply them and use all the hits.
	// The latter is particularly useful when directionality is added, as it can kill DR hits.
	fLimit_mPMT = true;

	fSTimePDFLimitsQueueNegative = -3.;				  //-5;//-3;//-1.5;//-3;//-10;
	fSTimePDFLimitsQueuePositive = 4.;				  // 10;//4;//3;//4;//500;
	fSTimePDFLimitsQueueNegative_fullTimeWindow = 0.; //-1.5;//-3;//-10;
	fSTimePDFLimitsQueuePositive_fullTimeWindow = 0.; // 3;//4;//500;
	fTimeWindowSizeFull = 1500;
	fAveraging = 20; // Number of points we used to average the vertex on.

	fMinimizeLimitsNegative = -700;
	fMinimizeLimitsPositive = 1000;
	
	fHitTimeLimitsNegative = -5; //-10;//-5;//-5    //-5
	fHitTimeLimitsPositive = 7;	 // 15;//7;//7		//7

	fSearchVtxStep = 300;		 // in cm, the step size for coarse grid search

	fSearchVtxTolerance = 60;	 // Number of candidate vertex that are kept after coarse grid search

	fIntegrationTimeWindow = 50; // in ns //? never Used ???
	DirTakeAll = false;

	// Compute light speed:
	float fCVacuum = 3e8 * 1e2 / 1e9; // speed of light, in centimeter per ns.
	float fNIndex = 1.373;			  // 1.385;//1.373;//refraction index of water
	fLightSpeed = fCVacuum / fNIndex;

	/*
	fTrueVtxPos.clear();
	fTrueVtxPosDouble.clear();
	fTrueVtxPos.resize(5,0.);
	fTrueVtxPosDouble.push_back(fTrueVtxPos);
	*/

	fPDFNorm_fullTimeWindow = 0.;
	
	std::map<std::string, double*> envVariables = 
	{
		{"SearchVtxStep", &fSearchVtxStep},
		{"SearchVtxTolerance", &fSearchVtxTolerance},
		{"HitTimeLimitsNegative", &fHitTimeLimitsNegative},
		{"HitTimeLimitsPositive", &fHitTimeLimitsPositive}
	};

	for (const auto& [envName, variable] : envVariables) 
	{
		const char* envValue = std::getenv(envName.c_str());
		if (envValue != nullptr) 
		{
			*variable = std::stod(envValue);
			std::cout << envName << " is set to " << *variable << std::endl;
		}
		else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
	}

}