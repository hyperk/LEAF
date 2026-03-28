/*********************************************************************************/
/**	analysis_fullrec.cpp							**/
/**	Author: Benjamin Quilain (benjamin.quilain@llr.in2p3.fr)		**/
/**	Date: April 10th 2025							**/
/**	Desc: Perform the vertex, direction and energy reconstruction	        **/
/*********************************************************************************/


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

#include "LEAF.hpp"	
#include "HKManager.hpp"	
//#define WITH_BONSAI
#ifdef WITH_BONSAI
#include "WCSimBonsai.hh"
#endif

//-----------------------------------------------------------------------------------------//

	int   eventId;					// event Id	
	int   triggerId;				// trigger Id (with event Id)

	double TotalQ;

	ROOT::Math::XYZVector particleDir(0,0,0);
	struct FitterAnalysis {
		double Wall;
		double Good;
		int n50[3];
		double dir[3][3];
		double dir_goodness[3];

		double dirKS[3];
	};

	void Normalize(double a[3]){
		double l;
		l=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
		for(int j=0; j<3; j++){
			a[j]/=l;
		}
	}

	//True information
	std::vector<int>   true_particleId; 		// True particle Id (PDGId)
	std::vector<double> true_energies;			// True energy (MeV)
	std::vector<ROOT::Math::XYZTVector> true_origins;		// True origin vertex (cm)
	std::vector<double> relativeAngle;
	std::vector<double> lf_relativeAngle;

	ROOT::Math::XYZVector lf_Dir;
	ROOT::Math::XYZVector lf_Dir_interp;
	double lf_Dir_res;
	WCSimRootGeom *geotree = 0; 

	//Raw hit
	std::vector<int>   rawhit_pmtId;  		// List of pmtId
	std::vector<int>   rawhit_Type;			// List of pmt Type
	std::vector<double> rawhit_T;			// List of hit times
	std::vector<double> rawhit_Q;			// List of hit charge'
	std::vector<bool>  rawhit_dark;  		// List of DarkNoise flag
	int rawhit_num;					// Number of rawhit
	int rawhit_num_noDN;				// Number of rawhit without DarkNoise


	std::vector<int>   rawhit_pmtId_tmp;  		// List of pmtId
	std::vector<int>   rawhit_Type_tmp;		// List of pmt Type
	std::vector<double> rawhit_T_tmp;		// List of hit times

	// In case of mPMT mixed with B&L we don't want to re apply the Digitization on the other type PMTs
	int fLastRawHit;


	//Digitized hit
	std::vector<int>   digithit_pmtId;		// List of pmtId
	std::vector<PMTType>   digithit_Type;		// List of pmt Type
	std::vector<double> digithit_T;			// List of hit times
	std::vector<double> digithit_Q;			// List of hit charge'
	std::vector<int> Charge_PMT;
	std::vector<bool>  digithit_dark;		// List of DarkNoise flag
	int digithit_num;					// Number of digit hit
	int digithit_num_noDN;				// Number of digit hit without DarkNoise

	//Bonsai output
	float bs_vertex[4];	   			// Bonsai reconstructed vertex (cm)
	float bs_good[3];				// Bonsai goodness
	float bs_energy;

	double fBSTime;
	double fLFTime;
	double lf_spatial_res;
	double lf_jitter;
	LEAFStructures::FitterOutput leaf_output;
	FitterAnalysis	leaf_output_ana;
	FitterAnalysis	bs_output_ana;

	int Hit_ID;
	int Hit_ID_50;
	int Hit_ID_200;
	int Hit_ID_400;

	int Hit_mPMT;
	int Hit_mPMT_50;
	int Hit_mPMT_200;
	int Hit_mPMT_400;

	int Hit_OD;
	int Hit_OD_50;
	int Hit_OD_200;
	int Hit_OD_400;

	int fHit     = 0;
	int fHit_50  = 0;
	int fHit_200 = 0;
	int fHit_400 = 0;
	double dWall =0;
	double lf_dWall =0;
	double lf_ToWall =0;
	WCSimRootGeom  * fGeometry;

	double R = 3242.96; 
	double h = 6701.41;

	// Input variable for Bonsai
	int bsCAB[2000];
	float bsT[2000];
	float bsQ[2000];
	int bsnhit[1]; 

//-----------------------------------------------------------------------------------------//


void SetCustomBranch(TTree* fPrimaryTree);
void SetCustomBranchInput(TTree* fPrimaryTree);
void SetGeoBranch(TTree* fGeoTree);

// Apply analysis on the event
bool AnalyseEvent(WCSimRootEvent * tEvent, PMTType iEventType, int i);
bool AnalyselfEvent(WCSimRootEvent * tEvent, PMTType iEventType);
TH1D *hRelativeAngle = nullptr;
TH1D *lf_hRelativeAngle = nullptr;
TH1D *hRelativeAngleCos = nullptr;
TH1D *lf_hRelativeAngleCos = nullptr;
TH1D *ToWall_Charge = nullptr;
TH2D *ToWall_Charge_2D = nullptr;
TH1D* Sum_HitQ = nullptr;
TH1D* Count_Hits = nullptr;

TH2D* hToWall_RelAngle_Charge = new TH2D("ToWall_RelAngle_Charge", "ToWall vs Rel. Angle and Charge", 100, 0, 9000, 360, 0, 180);
TProfile* pAngle_ToWall = new TProfile("pAngle_ToWall", "Angle vs ToWall", 360, 0, 180);

TH2D* hdWall_TotalCharge = new TH2D("DWall_TotalCharge", "dWall vs Total Charge", 150, 0, 3500, 50, 0, 400);
TH2D* hdWall_TotalCharge_corr = new TH2D("DWall_TotalChargeCorr", "dWall vs Corrected Total Charge", 150, 0, 3500, 50, 0, 400);
TH1D* hdWall_TotalCharge_hist = new TH1D("DWall_TotalCharge_hist", "dWall vs Total Charge (Hist)", 150, 0, 3500);
TH1D* Sum_HitQ_hdWall = new TH1D("Sum_HitQ", "Sum of HitQ per Wall bin", 150, 0, 3500);
TH1D* Count_Hits_hdWall = new TH1D("Count_Hits", "Count of Hits per Wall bin", 150, 0, 3500);



TH1D* hAngle_HitQ = new TH1D("Angle_HitQ", "Angle vs Hit Charge", 360, 0, 90);

void initializeHistograms() {
    hToWall_RelAngle_Charge->GetXaxis()->SetTitle("ToWall [cm]");
    hToWall_RelAngle_Charge->GetYaxis()->SetTitle("Relative Angle [degrees]");
    hToWall_RelAngle_Charge->GetZaxis()->SetTitle("Charge [p.e.]");

    hdWall_TotalCharge->GetXaxis()->SetTitle("dWall [cm]");
    hdWall_TotalCharge->GetYaxis()->SetTitle("Total Charge [p.e.]");

	hdWall_TotalCharge_corr->GetXaxis()->SetTitle("dWall [cm]");
    hdWall_TotalCharge_corr->GetYaxis()->SetTitle("Corrected Total Charge [p.e.]");

    hdWall_TotalCharge_hist->GetXaxis()->SetTitle("dWall [cm]");
    hdWall_TotalCharge_hist->GetYaxis()->SetTitle("Total Charge [p.e.]");

    hAngle_HitQ->GetXaxis()->SetTitle("Angle [degrees]");
    hAngle_HitQ->GetYaxis()->SetTitle("Hit Charge [p.e.]");
}

double calculateToWall(double R, double h, const ROOT::Math::XYZVector& position, const ROOT::Math::XYZVector& direction) {
    double x = position.X(), y = position.Y(), z = position.Z();
    double dx = direction.X(), dy = direction.Y(), dz = direction.Z();

    double A = dx * dx + dy * dy;
    double B = 2 * (x * dx + y * dy);
    double C = x * x + y * y - R * R;

    double discriminant = B * B - 4 * A * C;

    double t_cylinder = -1;
    if (discriminant >= 0) {
        double t1 = (-B - sqrt(discriminant)) / (2 * A);
        double t2 = (-B + sqrt(discriminant)) / (2 * A);

        if (t1 > 0) t_cylinder = t1;
        if (t2 > 0 && (t2 < t1 || t_cylinder < 0)) t_cylinder = t2;
    }

    //Calculate the t values for intersection with top and bottom
    double t_top = (h / 2 - z) / dz;
    double t_bottom = (-h / 2 - z) / dz;

    //Check if the intersections are on the top or bottom
    double x_top = x + t_top * dx;
    double y_top = y + t_top * dy;
    double x_bottom = x + t_bottom * dx;
    double y_bottom = y + t_bottom * dy;

    if (x_top * x_top + y_top * y_top > R * R) t_top = -1;       //Invalid if outside radius
    if (x_bottom * x_bottom + y_bottom * y_bottom > R * R) t_bottom = -1; //Invalid if outside radius

    //Find the smallest positive t that corresponds to a valid intersection (side or top/bottom)
    double t_wall = -1;
    if (t_cylinder > 0) t_wall = t_cylinder;
    if (t_top > 0 && (t_top < t_wall || t_wall < 0)) t_wall = t_top;
    if (t_bottom > 0 && (t_bottom < t_wall || t_wall < 0)) t_wall = t_bottom;

    return t_wall;
}

double calculateDWall(double R, double h, const ROOT::Math::XYZVector& position) {
    double x = position.X(), y = position.Y(), z = position.Z();
    double distanceToAxis = sqrt(x * x + y * y);
    double dSide = fabs(distanceToAxis - R);
    double dTop = fabs(z - h / 2);
    double dBottom = fabs(z + h / 2);
    return std::min({dSide, dTop, dBottom});
}






int main(int argc, char** argv){

	initializeHistograms();




	std::string sInputFile = "";
	std::string sOutputFile = "";

	int iNeededArgc = 5;	
	double dDarkNoise 	= 0.; //kHz
	double dDarkNoiseHybrid = 0.; //kHz
	iNeededArgc += 1;
	int StartEvent;
	int NumberEvents;
#ifdef mPMT	
	iNeededArgc += 1;
#endif

	if ( argc == iNeededArgc ) {
		int iArg = 1;
		sInputFile  	 =	argv[iArg];  iArg += 1;
		sOutputFile 	 =  	argv[iArg];  iArg += 1;
		dDarkNoise  	 = atof(argv[iArg]); iArg += 1;
#ifdef mPMT	
		dDarkNoiseHybrid = atof(argv[iArg]); iArg += 1;
#endif
	StartEvent = atoi(argv[iArg]); iArg += 1;
	NumberEvents = atoi(argv[iArg]); iArg += 1;

	}
	else {

		std::cout << "Synthax: " << argv[0] << " input output";
		std::cout << " DN_in_kHz_B&L";

#ifdef mPMT	
		std::cout << " DN_in_kHz_mPMT";
#endif
		std::cout << std::endl;
		return 0;
	}


	// Read WCSim output
	TFile* fInputFile = new TFile(sInputFile.c_str(),"READ");

	// Get TTrees
	TTree* fInputTree = (TTree*) fInputFile->Get("wcsimT");
	TTree *fInputGeoTree   = (TTree*) fInputFile->Get("wcsimGeoT");

	fGeometry = 0; 
	fInputGeoTree->SetBranchAddress("wcsimrootgeom",    &fGeometry);

	WCSimRootEvent * fIDevent  = new WCSimRootEvent();
	// Set Branche
	fInputTree->SetBranchAddress("wcsimrootevent" ,&fIDevent );
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent"   )->SetAutoDelete(kTRUE);

#ifdef mPMT
	WCSimRootEvent * fHybridevent  = new WCSimRootEvent();
	// Set Branche
	fInputTree->SetBranchAddress("wcsimrootevent2" ,&fHybridevent );
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent2"   )->SetAutoDelete(kTRUE);
#endif

#ifdef OD_ON
	WCSimRootEvent * fODevent  = new WCSimRootEvent();
	// Set Branche
	fInputTree->SetBranchAddress("wcsimrootevent_OD",&fODevent );
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);
#endif

	// Read Geo
	fInputGeoTree->GetEntry(0);


	// Create Output TTree
	TFile* fOutputFile    = new TFile(sOutputFile.c_str(),"RECREATE");
	fOutputFile->SetCompressionLevel(2);
	TTree* fGeoTree = new TTree("wcsimGeoT","Geometry TTree");	
	TTree* fPrimaryTree = new TTree("Reduced","Reduced TTree");

	// Set Branches
	SetGeoBranch(fGeoTree);
	SetCustomBranch(fPrimaryTree);

	fGeoTree->Fill();

	// Get PMT Number:
	int   nPMT_ID = fGeometry->GetWCNumPMT();
#ifdef OD_ON
	int   nPMT_OD = fGeometry->GetODWCNumPMT();
#endif
	int   nMultPMT = fGeometry->GetWCNumPMT(true);
	// TODO solve this issue!

	std::cout << " ID " << nPMT_ID << std::endl;
	std::cout << " mPMT " << nMultPMT << std::endl;

	HKManager::GetME()->SetGeometry(fGeometry,dDarkNoise * 1e3,dDarkNoiseHybrid * 1e3);

	// Initialize LEAF
	LEAF::GetME()->Initialize(	HKManager::GetME()->GetGeometry(),
								HKManager::GetME()->GetDarkNoise(),
								HKManager::GetME()->GetGeometryPMT_ID(),
								HKManager::GetME()->GetGeometryPMT_mPMT());
	LEAF::GetME()->SetNThread(); // Set number of Threads, default in the class is 12

	// Initialize HKAstroAnalysis
#ifdef WITH_HK_ASTROANALYSIS
	HKAstroAnalysis::GetME()->Initialize(HKManager::GetME()->GetGeometry());
#endif	

#ifdef WITH_BONSAI
	// Initialize Bonsai
	WCSimBonsai* fBonsai = new WCSimBonsai();
	fBonsai->Init(fGeometry);
#endif	
	// Read Input Tree
	int nPrimaryEvents = fInputTree->GetEntries();
	int iWrite         = 0;

	TStopwatch timer;
	timer.Reset();
	timer.Start();


	lf_hRelativeAngle = new TH1D("lf_ChargeProfile", "lf_ChargeProfile", 120, 0, 180);
	lf_hRelativeAngleCos = new TH1D("lf_ChargeProfileCos", "lf_ChargeProfileCos", 100, -1, 1);
	hRelativeAngle = new TH1D("ChargeProfile", "ChargeProfile", 120, 0, 180);
	hRelativeAngleCos = new TH1D("ChargeProfileCos", "ChargeProfileCos", 100, -1, 1);
	ToWall_Charge = new TH1D("ToWallCharge_hist", "ToWallCharge (Hist)", 1000, 0, 10000);
	ToWall_Charge_2D = new TH2D("ToWallCharge", "ToWallCharge", 1000, 0, 10000, 50, 0, 400);
	Sum_HitQ = new TH1D("Sum_HitQ", "Sum of HitQ per Wall bin", 1000, 0, 10000);
	Count_Hits = new TH1D("Count_Hits", "Count of Hits per Wall bin", 1000, 0, 10000);

	// Loop on Primary events	
	for(int i=StartEvent; i < StartEvent + NumberEvents; i++){
		// Reset Hit vector
		HKManager::GetME()->ResetHitInfo();

		//if ( i%1000==0 ) {
			timer.Stop();
			std::cout << "Event # = " << i << " / " << nPrimaryEvents << " ( " << timer.RealTime() << " )\n";
			timer.Reset();
			timer.Start();
		//}

		fInputTree->GetEntry(i); 

		// Initialize output variables
		eventId   		  = iWrite;
		triggerId 		  = 0;
		rawhit_num_noDN 	  = 0;

		true_particleId.clear();
		true_energies.clear();
		true_origins.clear();

		digithit_pmtId.clear();
		digithit_T.clear();
		digithit_Q.clear();
		Charge_PMT.clear();
		digithit_dark.clear();
		relativeAngle.clear();
		lf_relativeAngle.clear();
		rawhit_num 		  = 0;
		digithit_num		  = 0;

		bs_vertex	      [0] = -9999.;
		bs_vertex	      [1] = -9999.;
		bs_vertex	      [2] = -9999.;
		bs_vertex	      [3] = -9999.;
		bs_good		      [0] = -9999.;
		bs_good		      [1] = -9999.;
		bs_good		      [2] = -9999.;

		fBSTime      = -9999;
		fLFTime      = -9999;

		Hit_ID       = 0;
		Hit_ID_50    = 0;
		Hit_ID_200   = 0;
		Hit_ID_400   = 0;

		Hit_mPMT     = 0;
		Hit_mPMT_50  = 0;
		Hit_mPMT_200 = 0;
		Hit_mPMT_400 = 0;

		Hit_OD       = 0;
		Hit_OD_50    = 0;
		Hit_OD_200   = 0;
		Hit_OD_400   = 0;




#ifdef WITH_BONSAI
		// re-initialize Bonsai input
		for ( int iHit = 0; iHit < 2000; iHit++ ) {
			bsCAB[iHit] 	= 0;
			bsT  [iHit] 	= 0.;
			bsQ  [iHit] 	= 0.;
		}
		bsnhit[0] = 0;
#endif
		leaf_output.vtx = ROOT::Math::XYZTVector(0,0,0,0);
		leaf_output.vtx_nll	= -9999.;
		//leaf_output.InTime	= 0;
		leaf_output.energy = 0.;
		leaf_output.dir = ROOT::Math::XYZVector(0,0,0);
		//leaf_output.SNRList = std::vector<double>(); 

		leaf_output_ana.Wall	= 0;
		for ( int iType = 0; iType < 3; iType++ ) {
			leaf_output_ana.n50		[iType]		= 0;
			leaf_output_ana.dir		[iType][0]	= 0;
			leaf_output_ana.dir		[iType][1]	= 0;
			leaf_output_ana.dir		[iType][2]	= 0;
			leaf_output_ana.dir_goodness	[iType]	= 0;
			leaf_output_ana.dirKS		[iType]	= 0;
		}

		fLastRawHit = 0;

		/****************************************************************************************/
		/* ID events										*/
		/****************************************************************************************/

		//std::cout << " ID Event " << std::endl;
		fHit     = 0;
		fHit_50  = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		/*bool bID =*/ AnalyseEvent(fIDevent,PMTType::kID,i);



		Hit_ID     = fHit;
		Hit_ID_50  = fHit_50;
		Hit_ID_200 = fHit_200;
		Hit_ID_400 = fHit_400;

		//std::cout << " Hit: " << Hit_ID << " event " << i << std::endl;

		/****************************************************************************************/
		/* mPMT events										*/
		/****************************************************************************************/

		//std::cout << " mPMT Event " << std::endl;
		fHit     = 0;
		fHit_50  = 0;
		fHit_200 = 0;
		fHit_400 = 0;

#ifdef mPMT
		/*bool bmPMT =*/ AnalyseEvent(fHybridevent,PMTType::kmPMT,i);
#endif
		Hit_mPMT     = fHit;
		Hit_mPMT_50  = fHit_50;
		Hit_mPMT_200 = fHit_200;
		Hit_mPMT_400 = fHit_400;

		//std::cout << " Hit: " << Hit_mPMT << " event " << i << std::endl;

		/****************************************************************************************/
		/* OD events										*/
		/****************************************************************************************/
#ifdef OD_ON
		fHit     = 0;
		fHit_50  = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		// There should be a dedicated OD analyser, as many thing should be different than for ID
		// Doesn't exist yet
		/*bool bOD =*/ AnalyseODEvent(fODevent,PMTType::kOD);

		Hit_OD     = fHit;
		Hit_OD_50  = fHit_50;
		Hit_OD_200 = fHit_200;
		Hit_OD_400 = fHit_400;
#endif	
		/****************************************************************************************/
		/* Benjamin Fitter									*/
		/****************************************************************************************/

		//std::cout << " Start LEAF " << std::endl;
		TStopwatch timerLF;
		timerLF.Reset();
		timerLF.Start();

		LEAF::GetME()->LoadHitsCollection(HKManager::GetME()->GetHitCollection_ID(),HKManager::GetME()->GetHitCollection_mPMT());
        leaf_output = LEAF::GetME()->MakeSequentialFit();
        
		timerLF.Stop();

		fLFTime = timerLF.RealTime();

		//std::cout << " n50 " << leaf_output_ana.n50[0] << " " << leaf_output_ana.n50[1] << " " << leaf_output_ana.n50[2] << std::endl;
		std::cout << " LEAF took: " << timerLF.RealTime() << " for " << HKManager::GetME()->GetHitCollection_ID()->size() << " Hits"<<  std::endl;
		AnalyselfEvent(fIDevent,PMTType::kID);

#ifdef WITH_HK_ASTROANALYSIS
		HKAstroAnalysis::GetME()->SetVertex(leaf_output.Vtx);

		leaf_output_ana.Wall 		= HKAstroAnalysis::GetME()->ComputeDistanceFromWall();

		HKAstroAnalysis::GetME()->MakeAnalysis(HKManager::GetME()->GetHitCollection(),NormalPMT);

		leaf_output_ana.n50		[NormalPMT] 		= HKAstroAnalysis::GetME()->Getn50();
		leaf_output_ana.dirKS 		[NormalPMT]		= HKAstroAnalysis::GetME()->GetdirKS();
		leaf_output_ana.dir 		[NormalPMT][0]		= HKAstroAnalysis::GetME()->Getdir_Simple()[0];
		leaf_output_ana.dir 		[NormalPMT][1]		= HKAstroAnalysis::GetME()->Getdir_Simple()[1];
		leaf_output_ana.dir 		[NormalPMT][2]		= HKAstroAnalysis::GetME()->Getdir_Simple()[2];
		leaf_output_ana.dir_goodness	[NormalPMT]		= HKAstroAnalysis::GetME()->Getdir_Simple()[3];

		HKAstroAnalysis::GetME()->MakeAnalysis(HKManager::GetME()->GetHitCollection(),MiniPMT);

		leaf_output_ana.n50		[MiniPMT] 		= HKAstroAnalysis::GetME()->Getn50();
		leaf_output_ana.dirKS 		[MiniPMT]		= HKAstroAnalysis::GetME()->GetdirKS();
		leaf_output_ana.dir 		[MiniPMT][0]		= HKAstroAnalysis::GetME()->Getdir_Simple()[0];
		leaf_output_ana.dir 		[MiniPMT][1]		= HKAstroAnalysis::GetME()->Getdir_Simple()[1];
		leaf_output_ana.dir 		[MiniPMT][2]		= HKAstroAnalysis::GetME()->Getdir_Simple()[2];
		leaf_output_ana.dir_goodness	[MiniPMT]		= HKAstroAnalysis::GetME()->Getdir_Simple()[3];

		HKAstroAnalysis::GetME()->MakeAnalysis(HKManager::GetME()->GetHitCollection(),AllPMT);

		leaf_output_ana.n50		[AllPMT] 		= HKAstroAnalysis::GetME()->Getn50();
		leaf_output_ana.dirKS 		[AllPMT]		= HKAstroAnalysis::GetME()->GetdirKS();
		leaf_output_ana.dir 		[AllPMT][0]		= HKAstroAnalysis::GetME()->Getdir_Simple()[0];
		leaf_output_ana.dir 		[AllPMT][1]		= HKAstroAnalysis::GetME()->Getdir_Simple()[1];
		leaf_output_ana.dir 		[AllPMT][2]		= HKAstroAnalysis::GetME()->Getdir_Simple()[2];
		leaf_output_ana.dir_goodness	[AllPMT]		= HKAstroAnalysis::GetME()->Getdir_Simple()[3];

		leaf_output_ana.Good		= HKAstroAnalysis::GetME()->GoodnessBonsai();

#endif	

		if (true_origins.size() > 0 ) {
            ROOT::Math::XYZVector true_vertex(true_origins[0].X(), true_origins[0].Y(), true_origins[0].Z());
            ROOT::Math::XYZVector reconstructed_vertex(leaf_output.vtx.X(), leaf_output.vtx.Y(), leaf_output.vtx.Z());

            ROOT::Math::XYZVector  lf_vertex(reconstructed_vertex);
            ROOT::Math::XYZVector  lf_Dire(leaf_output.dir);
	    	lf_dWall = calculateDWall(R, h, lf_vertex);
		    lf_ToWall = calculateToWall(R, h, lf_vertex, lf_Dire);
		    lf_spatial_res = LEAFUtilities::GetDistance(true_vertex, reconstructed_vertex);
		    lf_jitter = abs(leaf_output.vtx.T()-true_origins[0].T());

		    TotalQ+=leaf_output.energy;
		    lf_Dir_interp = leaf_output.dir;
            double angle_lf = lf_Dir_interp.Dot(particleDir);


		lf_Dir_res = TMath::ACos(angle_lf)*180./TMath::Pi();


    } else {
        lf_spatial_res = -9999; // Assign invalid value if true vertex is missing
    }
#ifdef WITH_BONSAI
		/****************************************************************************************/
		/* Bonsai										*/
		/****************************************************************************************/

		if ( bsnhit[0] < 2000 && bsnhit[0] > 0 ) {

			TStopwatch timerBS;
			timerBS.Reset();
			timerBS.Start();

			// Some variable for Bonsai:
			float bsvertex[4];
			float bsresult[6];
			float bsgood[3];
			int bsnsel[1]; //nsel (SLE)

			// Fit with Bonsai
			//std::cout << "Bonsai hit: " << bsnhit[0] << std::endl;
			//for ( int iHit = 0; iHit < bsnhit[0]; iHit++ ) {
			//	std::cout << " hit " << iHit << " cab " << bsCAB[iHit] << " T " << bsT[iHit] << " Q " << bsQ[iHit] << std::endl;
			//}

			fBonsai->BonsaiFit( bsvertex, bsresult, bsgood, bsnsel, bsnhit, bsCAB, bsT, bsQ);

#ifdef WITH_HK_ASTROANALYSIS
			HKAstroAnalysis::GetME()->SetVertex(bs_vertex);

			bs_output_ana.Wall 		= HKAstroAnalysis::GetME()->ComputeDistanceFromWall();

			HKAstroAnalysis::GetME()->MakeAnalysis(HKManager::GetME()->GetHitCollection());

			bs_output_ana.n50		[NormalPMT] 		= HKAstroAnalysis::GetME()->Getn50();
			bs_output_ana.dirKS 		[NormalPMT]		= HKAstroAnalysis::GetME()->GetdirKS();
			bs_output_ana.dir 		[NormalPMT][0]		= HKAstroAnalysis::GetME()->Getdir_Simple()[0];
			bs_output_ana.dir 		[NormalPMT][1]		= HKAstroAnalysis::GetME()->Getdir_Simple()[1];
			bs_output_ana.dir 		[NormalPMT][2]		= HKAstroAnalysis::GetME()->Getdir_Simple()[2];
			bs_output_ana.dir_goodness	[NormalPMT]		= HKAstroAnalysis::GetME()->Getdir_Simple()[3];

			bs_output_ana.Good		= HKAstroAnalysis::GetME()->GoodnessBonsai();
#endif			

			bs_vertex[0] = bsvertex[0];
			bs_vertex[1] = bsvertex[1];
			bs_vertex[2] = bsvertex[2];
			bs_vertex[3] = bsvertex[3];

			bs_good[0]      = bsgood[0];
			bs_good[1]      = bsgood[1];
			bs_good[2]      = bsgood[2];

			//float diff = sqrt(pow(true_origin_X[0]- bsvertex[0],2.)+pow(true_origin_X[0]- bsvertex[0],2.)+pow(true_origin_X[0]- bsvertex[0],2.));

			//timerBS.Stop();
			//std::cout << " origin " << true_origin_X[0] << " " << true_origin_Y[0] << " " << true_origin_Z[0] << " diff: " << diff << std::endl;
			//std::cout << " BS took: " << timerBS.RealTime() << " goodness " << bsgood[0]<< " " << bsgood[1]<< " " <<  bsgood[2] << " bsvertex: " << bsvertex[0] << " " << bsvertex[1] << " " << bsvertex[2] << std::endl;


			fBSTime = timerBS.RealTime();
		}	
#endif
		/****************************************************************************************/
		/* Fill output tree									*/
		/****************************************************************************************/

		fPrimaryTree->Fill();
		iWrite += 1;  
	}
		std::cout << "Event # = " << nPrimaryEvents << " / " << nPrimaryEvents << std::endl;
// Convert histograms to graphs

fOutputFile->cd(); 





	lf_hRelativeAngle->Write(); 
	lf_hRelativeAngleCos->Write(); 
	hRelativeAngle->Write(); 
	hRelativeAngleCos->Write(); 
	ToWall_Charge->Write();
	hToWall_RelAngle_Charge->Write();
	pAngle_ToWall->Write();
	hdWall_TotalCharge->Write();
	hdWall_TotalCharge_corr->Write();
	hdWall_TotalCharge_hist->Write();
	ToWall_Charge_2D->Write();
	hAngle_HitQ->Write();

	fOutputFile->Write("",TObject::kOverwrite);

	// Cleaning
	delete fPrimaryTree;
	delete fOutputFile;

	//std::cout << "Total Charge: " << leaf_output.TotalCharge << std::endl;
	return 1;
//TotalQ 1000e10Mev: 207435 -- 100e10Mev: 20944.7 -- 1000e5MeV: 144119 -- 1000e3MeV: 45966.4
}

void SetCustomBranch(TTree* fPrimaryTree) {		

	fPrimaryTree->Branch("eventId",		&eventId,			"eventId/I");	
	fPrimaryTree->Branch("triggerId",		&triggerId,			"triggerId/I");

	fPrimaryTree->Branch("true_particleId",      &true_particleId);
	fPrimaryTree->Branch("true_origins",     	&true_origins);
	fPrimaryTree->Branch("true_energies",     	&true_energies);

	fPrimaryTree->Branch("rawhit_num", 		&rawhit_num,			"rawhit_num/I");
	fPrimaryTree->Branch("digithit_num", 		&digithit_num,			"digithit_num/I");

	fPrimaryTree->Branch("ID_hits", 		&Hit_ID,			"Hit_ID/I");	
	fPrimaryTree->Branch("ID_hits_50", 		&Hit_ID_50,			"Hit_ID_50/I");	
	fPrimaryTree->Branch("ID_hits_200", 		&Hit_ID_200,			"Hit_ID_200/I");	
	fPrimaryTree->Branch("ID_hits_400", 		&Hit_ID_400,			"Hit_ID_400/I");

	fPrimaryTree->Branch("mPMT_hits", 		&Hit_mPMT,			"Hit_mPMT/I");	
	fPrimaryTree->Branch("mPMT_hits_50", 		&Hit_mPMT_50,			"Hit_mPMT_50/I");	
	fPrimaryTree->Branch("mPMT_hits_200", 	&Hit_mPMT_200,			"Hit_mPMT_200/I");	
	fPrimaryTree->Branch("mPMT_hits_400", 	&Hit_mPMT_400,			"Hit_mPMT_400/I");

	fPrimaryTree->Branch("bs_vertex", 		bs_vertex,			"bs_vertex[4]/F");	
	fPrimaryTree->Branch("bs_good", 		bs_good,			"bs_good[3]/F");	
	fPrimaryTree->Branch("bs_ctime", 		&fBSTime,			"bs_ctime/F"); // Computation time

	fPrimaryTree->Branch("lf_spatial_res", &lf_spatial_res, "lf_spatial_res/D");
	fPrimaryTree->Branch("lf_jitter", &lf_jitter, "lf_jitter/D");
	fPrimaryTree->Branch("DWall", &dWall, "DWall/D");
	fPrimaryTree->Branch("lf_DWall", &lf_dWall, "DWall/D");
	fPrimaryTree->Branch("lf_ToWall", &lf_ToWall, "DWall/D");

	fPrimaryTree->Branch("reconstructed_energy", &leaf_output.energy, "reconstructed_energy/D");
	fPrimaryTree->Branch("TotalCharge", &leaf_output.total_charge, "TotalCharge/D");
	fPrimaryTree->Branch("ParticleDir", &particleDir,"particleDir[3]/D");
	fPrimaryTree->Branch("lf_Dir", &leaf_output.dir);
	fPrimaryTree->Branch("lf_Dir_res", &lf_Dir_res,"lf_Dir_res/D");
	//fPrimaryTree->Branch("SNR", &leaf_output.SNRList);

	fPrimaryTree->Branch("RelativeAngle", &relativeAngle);
	fPrimaryTree->Branch("lf_RelativeAngle", &lf_relativeAngle);

	fPrimaryTree->Branch("DigiHitQ", &digithit_Q);
	fPrimaryTree->Branch("Charge_PMT", &Charge_PMT);
	fPrimaryTree->Branch("DigiHitT", &digithit_T);




	fPrimaryTree->Branch("lf_vertex", 		&leaf_output.vtx);
	fPrimaryTree->Branch("lf_NLL", 		&leaf_output.vtx_nll,		"lf_NLL/D");
	//fPrimaryTree->Branch("lf_intime", 		&leaf_output.InTime,		"lf_intime/I");
	fPrimaryTree->Branch("lf_good", 		&leaf_output_ana.Good,		"lf_good/D");
	fPrimaryTree->Branch("lf_wall", 		&leaf_output_ana.Wall,		"lf_wall/D");
	fPrimaryTree->Branch("lf_n50", 		&leaf_output_ana.n50,		"lf_n50[3]/I");
	fPrimaryTree->Branch("lf_dir", 		&leaf_output_ana.dir,		"lf_dir[3][3]/D");
	fPrimaryTree->Branch("lf_dir_goodness", 	&leaf_output_ana.dir_goodness, "lf_dir_goodness[3]/D");
	fPrimaryTree->Branch("lf_dirKS", 		&leaf_output_ana.dirKS,	"lf_dirKS[3]/D");
	fPrimaryTree->Branch("lf_ctime", 		&fLFTime,			"lf_ctime/D"); // Computation time

}

void SetGeoBranch(TTree* fGeoTree) {	

	fGeoTree->Branch("wcsimrootgeom",		fGeometry);

}
bool AnalyselfEvent(WCSimRootEvent * tEvent, PMTType pmtType){

		int   nVertex			= 0; // Number of Vertex in event
		int   nTrack			= 0; // Number of Track in event
		int   nRawCherenkovHits	= 0; // Number of Raw Cherenkov hits
		//int   startingCherenkovHitID	= 0; // starting ID of Digitized Cherenkov hits. Usually starts at 0, but if there are two types of PMTs, it does not.
		int   nDigitizedCherenkovHits	= 0; // Number of Digitized Cherenkov hits
		//double 	fTriggerTime	= 0.;

		// Declare some useful variables
		WCSimRootTrigger * 	fRootTrigger;
		//TClonesArray *		fTimeArray;

		// Currently only one Trigger is used (to be check)
		//int iTrig = 0;
	for(int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++){

		triggerId    = iTrig;
		fRootTrigger = tEvent->GetTrigger(iTrig);

		// Grab the big arrays of times and parent IDs
		//fTimeArray   = fRootTrigger->GetCherenkovHitTimes();

		// Get number of vertex and tracks
#ifdef OLD_WCSIM
		nVertex		= 1;
#else
		nVertex		= fRootTrigger->GetNvtxs();
#endif

		nTrack			= fRootTrigger->GetNtrack();			
		nRawCherenkovHits	= fRootTrigger->GetNumTubesHit();

		if ( nTrack == 0 || nRawCherenkovHits ==  0 ) {
			// No track, no hit, nothing to do
			return false;
		}

		// Get number of hits
		//startingCherenkovHitID  = digithit_pmtId.size();
		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();	

		if (nDigitizedCherenkovHits == 0){
			continue;
		}
		//fTriggerTime = fRootTrigger->GetTriggerInfo()[2] - fRootTrigger->GetTriggerInfo()[1];
		ROOT::Math::XYZVector Truepos;
		// Loop on vertex
		for(int iVertex = 0; iVertex < nVertex; iVertex++){
#ifdef OLD_WCSIM
		    Truepos.SetXYZ(fRootTrigger->GetVtx(0),fRootTrigger->GetVtx(1),fRootTrigger->GetVtx(2));
#else
		    Truepos.SetXYZ(fRootTrigger->GetVtxs(iVertex,0),fRootTrigger->GetVtxs(iVertex,1),fRootTrigger->GetVtxs(iVertex,2));
#endif				
            if ( iVertex*2 > nTrack ) {
                std::cout << "ERROR: Vertex and Track number incompatible (nVertex: " << nVertex << " ; nTrack " << nTrack << ") " << std::endl;
                return 0;
            }

            // Beam info are registered in the Track
            TObject * Track = (fRootTrigger->GetTracks())->At(iVertex*2);
            WCSimRootTrack *wcTrack = dynamic_cast<WCSimRootTrack*>(Track);

            if (wcTrack->GetParentId() == 0 && wcTrack->GetIpnu()==11) {
                particleDir.SetXYZ(wcTrack->GetPdir(0),wcTrack->GetPdir(1),wcTrack->GetPdir(2));  
                particleDir = particleDir.Unit();
                }
    	    }

		    if ( nRawCherenkovHits < 1 ) return false;

			std::vector<double> times;

			// Get number of hit
			rawhit_num = nRawCherenkovHits;
			// Loop on Digitized Hit
            lf_Dir.SetXYZ(0,0,0);
			double AllQ = 0;
			for(int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
				ROOT::Math::XYZVector lf_particleRelativePMTpos;
				ROOT::Math::XYZVector particleRelativePMTpos;
				ROOT::Math::XYZVector lf_relativePMTpos;
				TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
				WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

				int pmtId      = wcDigitHit->GetTubeId();

				//double HitT     = wcDigitHit->GetT() + fTriggerTime;
				double HitQ     = wcDigitHit->GetQ();
				int peForTube = wcDigitHit->GetQ();



				WCSimRootPMT pmt;
				pmt = fGeometry->GetPMT(pmtId - 1, false);
				ROOT::Math::XYZVector PMTpos(pmt.GetPosition(0),pmt.GetPosition(1),pmt.GetPosition(2));
                particleRelativePMTpos = (PMTpos - Truepos);
                lf_relativePMTpos = (PMTpos - ROOT::Math::XYZVector(leaf_output.vtx.X(), leaf_output.vtx.Y(), leaf_output.vtx.Z()));
				ROOT::Math::XYZVector vDir(particleRelativePMTpos.Unit());
				ROOT::Math::XYZVector lfHitDir(lf_relativePMTpos.Unit());
                
                lf_Dir += (lfHitDir * HitQ);

				AllQ += HitQ;

				double hitlfRelativeAngle = vDir.Dot(leaf_output.dir);
                //std::cout << TMath::ACos(hitlfRelativeAngle)*180./TMath::Pi() << std::endl;

				lf_relativeAngle.push_back((TMath::ACos(hitlfRelativeAngle))*180./TMath::Pi());

				//std::cout << "Real Relative Angle : " << hitRelativeAngle*180./TMath::Pi() << std::endl;
				//std::cout << "Charge of hit : " << HitQ << std::endl;
				lf_hRelativeAngle->Fill((TMath::ACos(hitlfRelativeAngle))*180./TMath::Pi(), peForTube);
				lf_hRelativeAngleCos->Fill(hitlfRelativeAngle, peForTube);
			}
	}
	return true;

}
bool AnalyseEvent(WCSimRootEvent * tEvent, PMTType pmtType, int i) {

	int   nVertex			= 0; // Number of Vertex in event
	int   nTrack			= 0; // Number of Track in event
	int   nRawCherenkovHits	= 0; // Number of Raw Cherenkov hits
	int   startingCherenkovHitID	= 0; // starting ID of Digitized Cherenkov hits. Usually starts at 0, but if there are two types of PMTs, it does not.
	int   nDigitizedCherenkovHits	= 0; // Number of Digitized Cherenkov hits
	double 	fTriggerTime	= 0.;

	// Declare some useful variables
	WCSimRootTrigger * 	fRootTrigger;
	//TClonesArray *		fTimeArray;

	// Currently only one Trigger is used (to be check)
	//int iTrig = 0;
	for(int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++){

		triggerId    = iTrig;
		fRootTrigger = tEvent->GetTrigger(iTrig);

		// Grab the big arrays of times and parent IDs
		//fTimeArray   = fRootTrigger->GetCherenkovHitTimes();

		// Get number of vertex and tracks
#ifdef OLD_WCSIM
		nVertex		= 1;
#else
		nVertex		= fRootTrigger->GetNvtxs();
#endif

		nTrack			= fRootTrigger->GetNtrack();			
		nRawCherenkovHits	= fRootTrigger->GetNumTubesHit();

		if ( nTrack == 0 || nRawCherenkovHits ==  0 ) {
			// No track, no hit, nothing to do
			return false;
		}

		// Get number of hits

		startingCherenkovHitID  = digithit_pmtId.size();
		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();	

		if (nDigitizedCherenkovHits == 0){
		continue;
		}
		fTriggerTime            = fRootTrigger->GetTriggerInfo()[2] - fRootTrigger->GetTriggerInfo()[1];
		ROOT::Math::XYZVector Truepos;
		// Loop on vertex
		for(int iVertex = 0; iVertex < nVertex; iVertex++){
			// Beam info are registered in the Track
			TObject * Track = (fRootTrigger->GetTracks())->At(iVertex*2);
			WCSimRootTrack *wcTrack = dynamic_cast<WCSimRootTrack*>(Track);

#ifdef OLD_WCSIM
			Truepos.SetXYZ(fRootTrigger->GetVtx(0), fRootTrigger->GetVtx(1), fRootTrigger->GetVtx(2));
    		ROOT::Math::XYZTVector true_origin(fRootTrigger->GetVtx(0), fRootTrigger->GetVtx(1), fRootTrigger->GetVtx(2), wcTrack->GetTime());
            true_origins.push_back(true_origin);

#else
			Truepos.SetXYZ(fRootTrigger->GetVtxs(iVertex,0), fRootTrigger->GetVtxs(iVertex,1), fRootTrigger->GetVtxs(iVertex,2));
		    ROOT::Math::XYZTVector true_origin(fRootTrigger->GetVtxs(iVertex,0), fRootTrigger->GetVtxs(iVertex,1), fRootTrigger->GetVtxs(iVertex,2), wcTrack->GetTime());
		    true_origins.push_back(true_origin);
#endif				
			if ( iVertex*2 > nTrack ) {
                std::cout << "ERROR: Vertex and Track number incompatible (nVertex: " << nVertex << " ; nTrack " << nTrack << ") " << std::endl;
				return 0;
			}


            if (wcTrack->GetParentId() == 0 && wcTrack->GetIpnu()==11) {
                particleDir.SetXYZ(wcTrack->GetPdir(0),wcTrack->GetPdir(1),wcTrack->GetPdir(2));
                particleDir = particleDir.Unit();
            }

			leaf_output_ana.Wall = calculateToWall(R, h, Truepos, particleDir);
			dWall = calculateDWall(R, h, Truepos);
			std::cout << "true dir: " << particleDir.X() << "," << particleDir.Y() << "," << particleDir.Z() << std::endl;
			std::cout << "true pos: " << Truepos.X() << ","<< Truepos.Y() << ","<< Truepos.Z() << std::endl;
			//std::cout << "ToWall: " << leaf_output_ana.Wall << std::endl;
			true_particleId.push_back( wcTrack->GetIpnu() );
			true_energies.push_back(	   wcTrack->GetE()    );

		}

		// Send True Position to LEAF for check
		LEAF::GetME()->SetTrueVertexInfo(true_origins[0]);

		if ( nRawCherenkovHits < 1 ) return false;

		std::vector<double> times;

		// Get number of hit
		rawhit_num = nRawCherenkovHits;
		// Loop on Digitized Hit
        lf_Dir.SetXYZ(0,0,0);
		double AllQ = 0;



		for(int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
			ROOT::Math::XYZVector lf_particleRelativePMTpos;
			ROOT::Math::XYZVector particleRelativePMTpos;
			ROOT::Math::XYZVector lf_relativePMTpos;
			TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
			WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

			int pmtId      = wcDigitHit->GetTubeId();

			double HitT     = wcDigitHit->GetT() + fTriggerTime;
			double HitQ     = wcDigitHit->GetQ();
			int peForTube = wcDigitHit->GetQ();
			AllQ += HitQ;
			digithit_pmtId.push_back(pmtId);
			digithit_T.push_back(HitT);
			digithit_Q.push_back(HitQ);
			Charge_PMT.push_back(peForTube);

			Sum_HitQ->Fill(leaf_output_ana.Wall, HitQ);
			Count_Hits->Fill(leaf_output_ana.Wall, 1);

			Sum_HitQ_hdWall->Fill(dWall, HitQ);
			Count_Hits_hdWall->Fill(dWall, 1);

			WCSimRootPMT pmt;
			pmt = fGeometry->GetPMT(pmtId - 1, false);

			ROOT::Math::XYZVector PMTpos, PMTOrientation, vDir, lfHitDir;
			double normPMTOrientation = 0, NormvDir = 0, NormlfHitDir = 0;

            PMTpos.SetXYZ(pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2));
            
            PMTOrientation.SetXYZ(pmt.GetOrientation(0), pmt.GetOrientation(1), pmt.GetOrientation(2));
            PMTOrientation = PMTOrientation.Unit();

            vDir = (PMTpos - Truepos).Unit();
            lfHitDir = (PMTpos - ROOT::Math::XYZVector(leaf_output.vtx.X(),leaf_output.vtx.Y(),leaf_output.vtx.Z())).Unit();
            
			normPMTOrientation = sqrt(normPMTOrientation);
			NormvDir = sqrt(NormvDir);
			NormlfHitDir = sqrt(NormlfHitDir);

			double cosPMThitAngle = vDir.Dot(PMTOrientation);

			double PMThitAngle = acos(-cosPMThitAngle) * (180.0 / M_PI);
			//double cos2PMThitAngle = fabs(2 * cosPMThitAngle * cosPMThitAngle - 1);
			//double Ftheta = 0.205 + 0.524*cos2PMThitAngle + 0.390*(pow(cos2PMThitAngle,2)) - 0.132*(pow(cos2PMThitAngle,3));

			//std::cout << TMath::ACos(hitlfRelativeAngle)*180./TMath::Pi() << std::endl;
			double hitRelativeAngle = vDir.Dot(particleDir);

			relativeAngle.push_back((TMath::ACos(hitRelativeAngle))*180./TMath::Pi());
			digithit_dark.push_back(false);
			digithit_Type.push_back(pmtType);
			//std::cout << "Real Relative Angle : " << hitRelativeAngle*180./TMath::Pi() << std::endl;
			//std::cout << "Charge of hit : " << HitQ << std::endl;
			hRelativeAngle->Fill((TMath::ACos(hitRelativeAngle))*180./TMath::Pi(), peForTube);
			hRelativeAngleCos->Fill(hitRelativeAngle, peForTube);



			/*WallPMTAngle_Charge->Fill(leaf_output_ana.Wall, PMThitAngle*2, HitQ)*/

			hToWall_RelAngle_Charge->Fill(leaf_output_ana.Wall, (TMath::ACos(hitRelativeAngle))*180./TMath::Pi(), HitQ);
			pAngle_ToWall->Fill((TMath::ACos(hitRelativeAngle))*180./TMath::Pi(),leaf_output_ana.Wall);
			hAngle_HitQ->Fill(PMThitAngle, HitQ);

		}


		std::cout << "Hits: " << nDigitizedCherenkovHits << std::endl;
        lf_Dir /= AllQ;
        lf_Dir = lf_Dir.Unit();

		for (int bin = 1; bin <= ToWall_Charge->GetNbinsX(); ++bin) {
            double sum = Sum_HitQ->GetBinContent(bin);
            double count = Count_Hits->GetBinContent(bin);

            if (count > 0) {
                ToWall_Charge->SetBinContent(bin, sum / count); 
            } else {
                ToWall_Charge->SetBinContent(bin, 0); 
            }
		}
		ToWall_Charge_2D->Fill(leaf_output_ana.Wall, AllQ);

		for (int bin = 1; bin <= hdWall_TotalCharge_hist->GetNbinsX(); ++bin) {
            double sum = Sum_HitQ_hdWall->GetBinContent(bin);
            double count = Count_Hits_hdWall->GetBinContent(bin);

            if (count > 0) {
                hdWall_TotalCharge_hist->SetBinContent(bin, sum / count); 
            } else {
                hdWall_TotalCharge_hist->SetBinContent(bin, 0); 
            }
		}

		hdWall_TotalCharge->Fill(dWall, AllQ);
		hdWall_TotalCharge_corr->Fill(dWall, leaf_output.energy);

		/*double lf_rot_dir[3];
		lf_rot_dir[0]=lf_Dir[2];
		lf_rot_dir[1]=-lf_Dir[0];
		lf_rot_dir[2]=-lf_Dir[1];
		for(int j=0;j<3;j++){
			lf_Dir[j] = lf_rot_dir[j];
		}		
		lf_dir_res = sqrt(pow(lf_Dir[0] - particleDir[0], 2) + pow(lf_Dir[1] - particleDir[1], 2) + pow(lf_Dir[2] - particleDir[2], 2));*/

		int iIdx_BS = 0;

		nDigitizedCherenkovHits = digithit_pmtId.size();
		digithit_num = nDigitizedCherenkovHits;

		// Feed fitter
        if ( pmtType == PMTType::kID ) {
            for(int iDigitHit = startingCherenkovHitID; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
                times.push_back(digithit_T[iDigitHit]);

                //std::cout << " Add Hit with " << digithit_pmtId[iDigitHit] << " pmtType: " << pmtType << std::endl;
                if ( digithit_pmtId[iDigitHit] <= 0 ) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
                HKManager::GetME()->AddHit_ID(digithit_T[iDigitHit],digithit_Q[iDigitHit],digithit_pmtId[iDigitHit]);

                // Bonsai (do not store mPMT hits)
                if ( iIdx_BS < 2000 && digithit_T[iDigitHit] < 4000. && digithit_T[iDigitHit] > -4000. ) {
                    //std::cout << pmtId << " " << HitT << std::endl;
                    bsCAB[iIdx_BS] = digithit_pmtId[iDigitHit];
                    bsT  [iIdx_BS] = digithit_T[iDigitHit]; // shift BS time is needed if interaction time is < 0, this needs to be considered
                    bsQ  [iIdx_BS] = digithit_Q[iDigitHit];

                    iIdx_BS += 1;
                }
            }
        }
        else {
            for(int iDigitHit = startingCherenkovHitID; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
                times.push_back(digithit_T[iDigitHit]);
                //std::cout << " Add Hit with " << digithit_pmtId[iDigitHit] << " pmtType: " << pmtType << std::endl;
                if ( digithit_pmtId[iDigitHit] <= 0 ) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
                HKManager::GetME()->AddHit_mPMT(digithit_T[iDigitHit],digithit_Q[iDigitHit],digithit_pmtId[iDigitHit]);
            }
        }

		// Bonsai
		bsnhit[0]               = iIdx_BS;
		fHit			= digithit_pmtId.size();	

		// Compute hits
		std::vector<double> Hit_time_50;
		std::vector<double> Hit_time_200;
		std::vector<double> Hit_time_400;

		std::sort(times.begin(),times.end());
		for(int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){

			double HitT = times[iDigitHit];

			Hit_time_50.push_back(HitT);
			Hit_time_200.push_back(HitT);
			Hit_time_400.push_back(HitT);

			// Count hit in 50 ns window
			while ( HitT - Hit_time_50[0] > 50. )
				Hit_time_50.erase(Hit_time_50.begin());

			// Count hit in 200 ns window
			while ( HitT - Hit_time_200[0] > 200. )
				Hit_time_200.erase(Hit_time_200.begin());


			// Count hit in 400 ns window
			while ( HitT - Hit_time_400[0] > 400. )
				Hit_time_400.erase(Hit_time_400.begin());

			if ( (unsigned int) fHit_50 < Hit_time_50.size() ) fHit_50 = Hit_time_50.size();
			if ( (unsigned int) fHit_200 < Hit_time_200.size() ) fHit_200 = Hit_time_200.size();
			if ( (unsigned int) fHit_400 < Hit_time_400.size() ) fHit_400 = Hit_time_400.size();
		}

	}	


	return true;
}