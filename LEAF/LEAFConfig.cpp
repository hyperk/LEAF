#include "LEAFConfig.hpp"

std::mutex LEAFConfig::mtx = std::mutex();

double LEAFConfig::fTankRadius = 0.0; 
double LEAFConfig::fTankHeight = 0.0; 
double LEAFConfig::fTankHalfHeight = 0.0;

bool LEAFConfig::fStepByStep = false; // Step by step mode was not test and is not supported by multithreading
bool LEAFConfig::fUseDirectionality = false;
bool LEAFConfig::fHighEnergy = false;
bool LEAFConfig::fDoubleNLL = true;

// If true, we apply same cuts for timing as B&L PMT when searching vertex.
// Otherwise, we do not apply them and use all the hits.
// The latter is particularly useful when directionality is added, as it can kill DR hits.
bool LEAFConfig::fLimit_mPMT = true;

double LEAFConfig::fDirectionPDF_maxResidual = 15.0;
double LEAFConfig::fDirectionPDF_minResidual = -5.0;

double LEAFConfig::fSTimePDFLimitsQueueNegative = -3.;				  //-5;//-3;//-1.5;//-3;//-10;
double LEAFConfig::fSTimePDFLimitsQueuePositive = 4.;				  // 10;//4;//3;//4;//500;
double LEAFConfig::fSTimePDFLimitsQueueNegative_fullTimeWindow = 0.; //-1.5;//-3;//-10;
double LEAFConfig::fSTimePDFLimitsQueuePositive_fullTimeWindow = 0.; // 3;//4;//500;
double LEAFConfig::fTimeWindowSizeFull = 1500;
int LEAFConfig::fAveraging = 20; // Number of points we used to average the vertex on.

double LEAFConfig::fMinimizeLimitsNegative = -700;
double LEAFConfig::fMinimizeLimitsPositive = 1000;
	
double LEAFConfig::fHitTimeLimitsNegative = -5; //-10;//-5;//-5	//-5
double LEAFConfig::fHitTimeLimitsPositive = 7;	 // 15;//7;//7		//7

double LEAFConfig::fSearchVtxStep = 300;		 // in cm, the step size for coarse grid search

double LEAFConfig::fSearchVtxTolerance = 60;	 // Number of candidate vertex that are kept after coarse grid search

double LEAFConfig::fIntegrationTimeWindow = 50; // in ns //? never Used ???
bool LEAFConfig::fDirTakeAll = false;
double LEAFConfig::fThetaStep = 5 * TMath::DegToRad(); 
double LEAFConfig::fPhiStep = 5 * TMath::DegToRad(); 
int LEAFConfig::fDirTolerance = 50;
double LEAFConfig::fLightSpeed = C_WATER_CM_NS;
double LEAFConfig::fPDFNorm_fullTimeWindow = 0.;
std::vector<PMTType> LEAFConfig::fActivePMTTypes = {PMTType::kID};
//std::vector<PMTType> LEAFConfig::fActivePMTTypes = {PMTType::kID, PMTType::kmPMT};

void LEAFConfig::Initialize(double lTankRadius, double lTankHeight) {
	// Only write change here:
	
	fTankRadius = lTankRadius; 
	fTankHeight = lTankHeight; 
	fTankHalfHeight = fTankHeight / 2.;

	//*** These maps describe what will be loaded from the config file given */

	std::map<std::string, double*> envVariables = {
		{"SearchVtxStep", &fSearchVtxStep},
		{"SearchVtxTolerance", &fSearchVtxTolerance},
		{"HitTimeLimitsNegative", &fHitTimeLimitsNegative},
		{"HitTimeLimitsPositive", &fHitTimeLimitsPositive},
		{"MinimizeLimitsNegative", &fMinimizeLimitsNegative},
		{"MinimizeLimitsPositive", &fMinimizeLimitsPositive},
		{"IntegrationTimeWindow", &fIntegrationTimeWindow},
		{"TimeWindowSizeFull", &fTimeWindowSizeFull},
		{"STimePDFLimitsQueueNegative", &fSTimePDFLimitsQueueNegative},
		{"STimePDFLimitsQueuePositive", &fSTimePDFLimitsQueuePositive},
		{"STimePDFLimitsQueueNegative_fullTimeWindow", &fSTimePDFLimitsQueueNegative_fullTimeWindow},
		{"STimePDFLimitsQueuePositive_fullTimeWindow", &fSTimePDFLimitsQueuePositive_fullTimeWindow},
		{"DirectionPDF_maxResidual", &fDirectionPDF_maxResidual},
		{"DirectionPDF_minResidual", &fDirectionPDF_minResidual},
		{"ThetaStep", &fThetaStep},
		{"PhiStep", &fPhiStep}
	};

	std::map<std::string, int*> envIntVariables = {
		{"Averaging", &fAveraging},
		{"DirTolerance", &fDirTolerance}
	};

	std::map<std::string, bool*> envBoolVariables = {
		{"Limit_mPMT", &fLimit_mPMT},
		{"UseDirectionality", &fUseDirectionality},
		{"StepByStep", &fStepByStep},
		{"HighEnergy", &fHighEnergy},
		{"DoubleNLL", &fDoubleNLL},
		{"DirTakeAll", &fDirTakeAll},
	};

	for (const auto& [envName, variable] : envVariables) {
		const char* envValue = std::getenv(envName.c_str());
		if (envValue != nullptr) {
			*variable = std::stod(envValue);
			std::cout << envName << " is set to " << *variable << std::endl;
		}
		else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
	}

	for (const auto& [envName, variable] : envIntVariables) {
		const char* envValue = std::getenv(envName.c_str());
		if (envValue != nullptr) {
			*variable = std::stoi(envValue);
			std::cout << envName << " is set to " << *variable << std::endl;
		}
		else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
	}

	for (const auto& [envName, variable] : envBoolVariables) {
		const char* envValue = std::getenv(envName.c_str());
		if (envValue != nullptr) {
			std::string valueStr(envValue);
			if (valueStr == "true" || valueStr == "1") {
				*variable = true;
				std::cout << envName << " is set to true" << std::endl;
			} 
			else if (valueStr == "false" || valueStr == "0") {
				*variable = false;
				std::cout << envName << " is set to false" << std::endl;
			} 
			else std::cout << "Environment variable " << envName << " has an invalid value. Using default value : "<< *variable << std::endl;
		}
		else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
	}

}