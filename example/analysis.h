
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TStopwatch.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimEnumerations.hh"

#include "LEAF.hh"	
#include "LeafDefinitions.hh"
#include "HKManager.hh"	

#define OLD_WCSIM // To be used if WCSim version is older than 1.8 (i.e. without multi vertex)
// #define mPMT // To be used if you are using mPMT

#define ID_EVENT 1
#define OD_EVENT 2
// #define mPMT_EVENT	0
#define UNDEFINED_EVENT 0

//-----------------------------------------------------------------------------------------//

WCSimRootGeom *geotree = 0; 
WCSimRootGeom *fLeafGeometry;

//* Event True Properties
int eventId;   // event Id
int nTrigger; // trigger Id (with event Id)
int usedTriggerId; // trigger Id used for the fit
double fLfTriggerTime;
int bestTrigger;
double rawTriggerTime;
std::vector<int> true_particleId;  // True particle Id (PDGId)
std::vector<double> true_energy;   // True energy (MeV)
std::vector<double> true_origin_X; // True origin vertex (cm)
std::vector<double> true_origin_Y; // True origin vertex (cm)
std::vector<double> true_origin_Z; // True origin vertex (cm)
std::vector<double> true_origin_T; // True origin vertex (ns)
std::vector<double> trueDir;
std::vector<double> trueVertex;
double dWall =0;
double toWall =0;
int digithit_num;				 // Number of digit hit
int digithit_num_noDN;			 // Number of digit hit without DarkNoise
int hit_5_15ns; // number of hits with residual time between -5 and 15ns 
int hit_50ns; // number of hits with residual time between -50 and 50ns
int fLastRawHit; // In case of mPMT mixed with B&L we don't want to re apply the Digitization on the other type PMTs
int startingCherenkovHitID = 0;	 // starting ID of Digitized Cherenkov hits. Usually starts at 0, but if there are two types of PMTs, it does not.

//* Properties Estimated by Leaf
std::vector<double> leaf_Vertex; // the estimated vertex
std::vector<double> leaf_Dir; // the estimated direction
std::vector<double> leaf_MyDir; // estimated direction with one prior that is the the quick direction
std::vector<double> leaf_QuickDir; // quick computed directions from to_hit vectors

//* Resulting properties
double lf_spatial_res;
double lf_time_res;
double lf_Dir_res;
double lf_Quick_Dir_res;
double lf_MyDir_res;
double lf_energy_res;
double lf_dWall =0;
double lf_ToWall =0;

//* Raw hit infos
std::vector<int> rawhit_pmtId; // List of pmtId
std::vector<int> rawhit_Type;  // List of pmt Type
std::vector<double> rawhit_T;  // List of hit times
std::vector<double> rawhit_Q;  // List of hit charge'
std::vector<bool> rawhit_dark; // List of DarkNoise flag
int rawhit_num;				   // Number of rawhit
int rawhit_num_noDN;		   // Number of rawhit without DarkNoise
std::vector<int> rawhit_pmtId_tmp; // List of pmtId
std::vector<int> rawhit_Type_tmp;  // List of pmt Type
std::vector<double> rawhit_T_tmp;  // List of hit times

//* Digitized hit infos
std::vector<int> digithit_pmtId; // List of pmtId
std::vector<int> digithit_Type;	 // List of pmt Type
std::vector<double> digithit_T;	 // List of hit times
std::vector<double> correctedDigithit_T;
std::vector<bool> hit_is_DR;
std::vector<double> digithit_Q;	 // List of hit charge'
std::vector<double> digithit_Angle;	 // List of hit charge'
std::vector<double> digithit_Angle_NLL; // todo delete this
std::vector<double> digithit_NormAngle;	 // List of hit charge'
std::vector<double> digithit_NeighborsDist;  // List of hit charge'
std::vector<double> relativeAngle;
std::vector<double> lf_relativeAngle;
std::vector<int> Charge_PMT;
std::vector<double> hit_residual; // List of residual times

//* Computation Times
FitterOutputProps fOutputProps;
double fLFTime;

//* Other hits infos
int Hit_ID;
int Hit_ID_20;
int Hit_ID_50;
int Hit_ID_200;
int Hit_ID_400;

int Hit_mPMT;
int Hit_mPMT_20;
int Hit_mPMT_50;
int Hit_mPMT_200;
int Hit_mPMT_400;

int Hit_OD;
int Hit_OD_50;
int Hit_OD_200;
int Hit_OD_400;

int fHit = 0;
int fHit_20 = 0;
int fHit_50 = 0;
int fHit_200 = 0;
int fHit_400 = 0;

//* Constants
double R = 3242.96; 
double h = 6701.41;

//* Leaf parameters
double maxHitAngle = 190.0;
double N_Neighbors = 5;
double maxDistanceToNeighbors = 8000.0;

//* Input variables for Bonsai
int bsCAB[2000];
float bsT[2000];
float bsQ[2000];
int bsnhit[1];


//-----------------------------------------------------------------------------------------//

struct arguments
{
    char * inputFile  = NULL;
    char * outputFile = NULL;
	double darkNoise  = 4.2;	// Dark noise frequency in Hz
	double darkNoiseH = 0.;		// Dark noise frequency in Hz when using the hybrid geometry
	double timeshift  = 0.;		// Shift hit time
    int startEvent    = 0;		// First event to analyze
    int endEvent      = 0;		// Last event to analyze
	bool hybrid       = false;
    bool verbose      = false;
};

struct FitterAnalysis 
{
	double Wall;
	double Good;
	int n50[3];
	double dir[3][3];
	double dir_goodness[3];
	double dirKS[3];
};


//-----------------------------------------------------------------------------------------//

void SetCustomBranch(TTree *fPrimaryTree, FitterOutput leaf_output);
void SetCustomBranchInput(TTree *fPrimaryTree);
void SetGeoBranch(TTree *fGeoTree);
bool AnalyseEvent(WCSimRootEvent *tEvent, int iEventType);
bool PostLeafAnalysis(WCSimRootEvent * tEvent, int iEventType, FitterOutput leaf_output);
arguments FetchInput(int argc, char* argv[]);



