/*********************************************************************************/
/**	HKAstroAnalysis.hh							   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: June 17th 2020							   **/
/**	Desc: Manage fitter and PMT hits for HK analysis			   **/
/*********************************************************************************/


#ifndef HKAstroAnalysis_hh
#define HKAstroAnalysis_hh

//#include <algorithm>
//#include <ctime>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>

#include "HKManager.hh"

//ROOT Headers
#include "TMath.h"


#define CNS2CM 			21.58333

class HKAstroAnalysis/* : public TObject */{

	public:
		HKAstroAnalysis();
		~HKAstroAnalysis();
		static HKAstroAnalysis* GetME();
				
			
		void MakeAnalysis(int iType=NormalPMT);	
		void ResetAnalysis();
		
		void SetVertex(std::vector<double> lVtx);	
		void SetVertex(float* lVtx);	
		
		unsigned int		Getn20() 		{	return fInTime20.size(); 	}
		unsigned int 		Getn30() 		{	return fInTime30.size(); 	}
		unsigned int 		Getn50() 		{	return fInTime50.size(); 	}
		
		double 		GetdirKS() 		{	return fVtxReco_dirKS; 	}
		std::vector<double> 	Getdir_Simple() 	{	return fVtxReco_dirSimple; 	}
		std::vector<double> 	Getdir()	 	{	return fVtxReco_dir; 		} // Not set currently
		
		double 		Getdir_goodness()	{ 	return fVtxReco_dirSimple[3];	}
				
		
		// Get distance between final vertex and wall:
		double 		ComputeDistanceFromWall();
		
		// Return Bonsai-like goodness of fit
		double 		GoodnessBonsai();
		
		
		
	private:
		// WCSim objects:
		HKManager* myManager;
		static HKAstroAnalysis* myAna;

		// Input:
		std::vector<double> fRecoVtx;
		
		double fRecoVtx_X; 
		double fRecoVtx_Y; 
		double fRecoVtx_Z; 
		double fRecoVtx_T; 
		double fRecoVtx_R2; 
		double fRecoVtx_R;
		double fRecoVtx_Good;
		
		// Analysis:
		// Wall
		double fRecoVtx_Wall;
		
		// Directions
		std::vector<double> fVtxReco_dirSimple;
		std::vector<double> fVtxReco_dir;
		
		// Energy
		double fRecoVtx_E;
		double fRecoVtx_E_Simple;
		
		// DirKS
		double fVtxReco_dirKS;
		double fVtxReco_good;
		
		
		struct PMTHitExt {
			double T;
			double Q;
			int PMT;
			
			double Dist;
			double ToF;
			
			double NormX;
			double NormY;
			double NormZ;
			
			double theta;
			double phi;
		};
		
		struct SortingToF { 
			bool operator()(PMTHitExt const &a, PMTHitExt const &b) const { 
				return a.ToF < b.ToF;
			}
		};
			
		// InTime hits
		std::vector< PMTHitExt > fInTime20;
		std::vector< PMTHitExt > fInTime30;
		std::vector< PMTHitExt > fInTime50;
		std::vector< PMTHitExt > fHitExtInfo;
		
	private:
		// Analysis function
	
		// Load some constants needed to compute the direction
		void LoadDirConstant();						

		// Geometric function
		double GroupPMTs(int pmtType, int pmt_number);
		
		// Get number of in time hits in twindow from tof subtracted t-distribution 
		void ComputeInTimeHit(int iType=NormalPMT);
		
		// Return the detection efficiency for a given incident angle
		double GetPMTAngleEfficiency(double dCosPM);
		
		// Return a vector from the direction and the angles
		std::vector<double> TransformInVector(std::vector<double> lDir, double dPhi, double dCosTheta );
		
		// Return dirKS
		double GetDirKS(std::vector<double> lDir);
		
		// Return ?
		double GetDir_lkflf(double dCos, int iE, bool bSimple);
		
		// Return the direction
		// if Simple = false, energy is ignored
		// if Simple = true, energy is needed
		std::vector<double> GetDirection(double dEnergy, std::vector<PMTHitExt> tInTime, bool bSimple, double dCutOff);
		
		// Constant for direction fit
		double fEP10[21];
		std::vector< std::vector<double> > fHitPat;
		
		// Constant for EffHit computation
		double fEffHitLambda[12];
		int fNearPMT;
		double fOccupancyTable[10][10];
		
};

#endif

