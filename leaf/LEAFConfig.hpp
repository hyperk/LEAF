#ifndef LEAFCONFIG_HPP
#define LEAFCONFIG_HPP

#include <iostream>
#include <map>
#include <vector>
#include <TMath.h>
#include "Enums/PMTType.hpp"

#define VERBOSE 		0

#undef VERBOSE_VTX // In SearchVertex
#undef VERBOSE_NLL // In Likelihood Computation

#define VTX_X			0
#define VTX_Y			1
#define VTX_Z			2
#define VTX_T			3

#define N_THREAD		2 // Default for sukap
#define C_WATER_CM_NS   21.849964 // light celerity in water, in centimeter per ns.

namespace LEAFConfig {

    //* Geometry & Constants
    double fTankRadius;
    double fTankHeight;
    double fTankHalfHeight;
    double fLightSpeed = C_WATER_CM_NS;

    //* PDFs
    double fSTimePDFLimitsQueueNegative = -3.;
    double fSTimePDFLimitsQueuePositive = 4.;
    double fSTimePDFLimitsQueueNegative_fullTimeWindow = 0.;
    double fSTimePDFLimitsQueuePositive_fullTimeWindow = 0.;
    double fPDFNorm_fullTimeWindow = 0.;

    double fDirectionPDF_maxResidual = 15.0;
    double fDirectionPDF_minResidual = -5.0;

    //* Overall Fit Parameters
    double fTimeWindowSizeFull = 1500.;
    double fIntegrationTimeWindow = 50.;
    bool   fStepByStep = false;
    bool   fUseDirectionality = false;
    bool   fDoubleNLL = false;
    bool   fLimit_mPMT = true;
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
    bool fDirTakeAll = false;
    double fThetaStep = 5 * TMath::DegToRad(); 
    double fPhiStep = 5 * TMath::DegToRad(); 
    int fDirTolerance = 50;

    //! List of active PMT types
    std::vector<PMTType> fActivePMTTypes; 

    void Initialize(double lTankRadius, double lTankHeight) {
        fTankRadius = lTankRadius; 
        fTankHeight = lTankHeight; 
        fTankHalfHeight = fTankHeight / 2.;

        fStepByStep = false; // Step by step mode was not test and is not supported by multithreading
        fUseDirectionality = false;
        fHighEnergy = false;
        fDoubleNLL = true;

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
        fDirTakeAll = false;
        fThetaStep = 5 * TMath::DegToRad(); 
        fPhiStep = 5 * TMath::DegToRad(); 
        fDirTolerance = 50;

        //* Compute light speed:
        fLightSpeed = C_WATER_CM_NS;

        fPDFNorm_fullTimeWindow = 0.;

        fActivePMTTypes = {PMTType::kID};
        //fActivePMTTypes = {PMTType::kID, PMTType::kmPMT};
        
        //*** These maps describe what will be loaded from the config file given */

        std::map<std::string, double*> envVariables = 
        {
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

        std::map<std::string, int*> envIntVariables = 
        {
            {"Averaging", &fAveraging},
            {"DirTolerance", &fDirTolerance}
        };

        std::map<std::string, bool*> envBoolVariables = 
        {
            {"Limit_mPMT", &fLimit_mPMT},
            {"UseDirectionality", &fUseDirectionality},
            {"StepByStep", &fStepByStep},
            {"HighEnergy", &fHighEnergy},
            {"DoubleNLL", &fDoubleNLL},
            {"DirTakeAll", &fDirTakeAll},
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

        for (const auto& [envName, variable] : envIntVariables) 
        {
            const char* envValue = std::getenv(envName.c_str());
            if (envValue != nullptr) 
            {
                *variable = std::stoi(envValue);
                std::cout << envName << " is set to " << *variable << std::endl;
            }
            else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
        }

        for (const auto& [envName, variable] : envBoolVariables) 
        {
            const char* envValue = std::getenv(envName.c_str());
            if (envValue != nullptr) 
            {
                std::string valueStr(envValue);
                if (valueStr == "true" || valueStr == "1") 
                {
                    *variable = true;
                    std::cout << envName << " is set to true" << std::endl;
                } 
                else if (valueStr == "false" || valueStr == "0") 
                {
                    *variable = false;
                    std::cout << envName << " is set to false" << std::endl;
                } 
                else std::cout << "Environment variable " << envName << " has an invalid value. Using default value : "<< *variable << std::endl;
            }
            else std::cout << "Environment variable " << envName << " is not set. Using default value : "<< *variable << std::endl;
        }

    }
}; // namespace LEAFConfig

#endif