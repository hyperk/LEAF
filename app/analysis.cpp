/*********************************************************************************/
/**	analysis.cpp								**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		**/
/**	Date: Febuary 10th 2020							**/
/**	Desc: Example application code for Benjamin's Low-E Fitter for Hyper-K	**/
/*********************************************************************************/
/**	Modified: Benjamin Quilain (benjamin.quilain@llr.in2p3.fr)		**/
/**	Date: May 14th 2020							**/
/**	Desc: Modified the output format and allow run w/o BONSAI	        **/
/*********************************************************************************/


#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>

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

#include "BQFitter.hh"
//#define BONSAI // To be used if we wish to also use BONSAI as a fitter. Requires BONSAI libraries (see GNUMakefile)
#define OLD_WCSIM // To be used if WCSim version is older than 1.8 (i.e. without multi vertex)
#define mPMT // To be used if you are using mPMT
#define V2BQDIR // To be validated before pushing
#ifdef BONSAI
#include "WCSimBonsai.hh"
#endif

#define ID_EVENT 	1
#define OD_EVENT	2
#define mPMT_EVENT	3
#define UNDEFINED_EVENT	0

//-----------------------------------------------------------------------------------------//

	int   eventId;					// event Id	
	int   triggerId;				// trigger Id (with event Id)

	struct Hit {
		int pmtId;
		float T;
		float Q;
		bool dark;
		int Type;
	};
	struct SortingT { 
		bool operator()(Hit const &a, Hit const &b) const { 
			return a.T < b.T;
		}
	};

	//True information
	std::vector<int>   true_particleId; 		// True particle Id (PDGId)
	std::vector<float> true_energy;			// True energy (MeV)
        std::vector<float> true_vertex_X;		// True origin vertex (cm)
        std::vector<float> true_vertex_Y;		// True origin vertex (cm)
        std::vector<float> true_vertex_Z;		// True origin vertex (cm)
        std::vector<float> true_vertex_T;		// True origin vertex time (ns)
        std::vector<float> true_dir_X;		        // True direction (unit)
        std::vector<float> true_dir_Y;		        // True direction (unit)
        std::vector<float> true_dir_Z;		        // True direction (unit)
		
	//Raw hit
	std::vector<int>   rawhit_pmtId;  		// List of pmtId
	std::vector<int>   rawhit_Type;			// List of pmt Type
	std::vector<float> rawhit_T;			// List of hit times
	std::vector<float> rawhit_Q;			// List of hit charge'
	std::vector<bool>  rawhit_dark;  		// List of DarkNoise flag
	int rawhit_num;					// Number of rawhit
	int rawhit_num_noDN;				// Number of rawhit without DarkNoise
	
	
	std::vector<int>   rawhit_pmtId_tmp;  		// List of pmtId
	std::vector<int>   rawhit_Type_tmp;		// List of pmt Type
	std::vector<float> rawhit_T_tmp;		// List of hit times
	
	// In case of mPMT mixed with B&L we don't want to re apply the Digitization on the other type PMTs
	int fLastRawHit;
		
	
	//Digitized hit
	std::vector<int>   digithit_pmtId;		// List of pmtId
	std::vector<int>   digithit_Type;		// List of pmt Type
	std::vector<float> digithit_T;			// List of hit times
	std::vector<float> digithit_Q;			// List of hit charge'
	std::vector<bool>  digithit_dark;		// List of DarkNoise flag
	int digithit_num;					// Number of digit hit
	int digithit_num_noDN;				// Number of digit hit without DarkNoise
	
	//Bonsai output
	float bs_vertex[4];	   			// Bonsai reconstructed vertex (cm)
	float bs_good[3];				// Bonsai goodness
	float bs_energy;
	
	float fBSTime;
	float fBQTime;

	BQFitter::FitterOutput leaf_output;
	
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
	
	WCSimRootGeom  * fGeometry;
				
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
bool AnalyseEvent(WCSimRootEvent * tEvent, int iEventType);

int main(int argc, char** argv){

	std::string sInputFile = "";
	std::string sOutputFile = "";

	int iNeededArgc = 5;
	int iStartEvent = 0;//ID of the first event to fit
	int iEndEvent   = 0;//ID of the final event to fit
	double dDarkNoise 	= 0.; //kHz
	double dDarkNoiseHybrid = 0.; //kHz
	iNeededArgc += 1;
#ifdef mPMT	
	iNeededArgc += 1;
#endif
	
	if ( argc == iNeededArgc ) {
		int iArg = 1;
		sInputFile  	 =	argv[iArg];  iArg += 1;
		sOutputFile 	 =  	argv[iArg];  iArg += 1;
		iStartEvent      = atoi(argv[iArg]);  iArg += 1;
		iEndEvent        = atoi(argv[iArg]);  iArg += 1;
		dDarkNoise  	 = atof(argv[iArg]); iArg += 1;
#ifdef mPMT	
		dDarkNoiseHybrid = atof(argv[iArg]); iArg += 1;
#endif
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
	fInputTree->SetBranchAddress("wcsimrootevent"   ,&fIDevent );
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent"   )->SetAutoDelete(kTRUE);
	
	WCSimRootEvent * fHybridevent  = new WCSimRootEvent();
#ifdef mPMT
	// Set Branche
	fInputTree->SetBranchAddress("wcsimrootevent2"   ,&fHybridevent );
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
	
	std::cout << " ID " << nPMT_ID << std::endl;
	std::cout << " mPMT " << nMultPMT << std::endl;

	// Threading is not supported on Sukap jobs
	// On sukap machine 10 cores is the maximum allowed (12 thread allow to use about 10 cores at maximum)
	BQFitter::GetME()->SetNThread(12);
	BQFitter::GetME()->SetGeometry(fGeometry,dDarkNoise * 1e3,dDarkNoiseHybrid * 1e3);

#ifdef BONSAI
	// Initialize Bonsai
	WCSimBonsai* fBonsai = new WCSimBonsai();
	fBonsai->Init(fGeometry);
#endif
	// Read Input Tree
	int nPrimaryEvents = std::min((Long64_t) iEndEvent,fInputTree->GetEntries());
	//int nPrimaryEvents = TMath::Min((Long64) iEndEvent,fInputTree->GetEntries());
	int iWrite         = 0;
	
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	
	// Loop on Primary events	
	for(int i=iStartEvent; i < nPrimaryEvents; i++){ 
	
		// Reset Hit vector
		BQFitter::GetME()->ResetHitInfo();
		
		if ( i%1000==0 ) {
			timer.Stop();
			std::cout << "Event # = " << i << " / " << (nPrimaryEvents - iStartEvent) << " ( " << timer.RealTime() << " )\n";
			timer.Reset();
			timer.Start();
		}
	
		fInputTree->GetEvent(i); 
		
		// Initialize output variables
		eventId   		  = iWrite;
		triggerId 		  = 0;
		rawhit_num_noDN 	  = 0;
		
		true_particleId.clear();
		true_energy.clear();
		true_vertex_X.clear();
		true_vertex_Y.clear();
		true_vertex_Z.clear();
		true_vertex_T.clear();
		true_dir_X.clear();
		true_dir_Y.clear();
		true_dir_Z.clear();
				
		digithit_pmtId.clear();
		digithit_T.clear();
		digithit_Q.clear();
		digithit_dark.clear();
		
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
		fBQTime      = -9999;
				
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


		// re-initialize Bonsai input
		for ( int iHit = 0; iHit < 2000; iHit++ ) {
			bsCAB[iHit] 	= 0;
			bsT  [iHit] 	= 0.;
			bsQ  [iHit] 	= 0.;
		}
		bsnhit[0] = 0;

		leaf_output.Vtx[0]	= 0.;
		leaf_output.Vtx[1]	= 0.;
		leaf_output.Vtx[2]	= 0.;
		leaf_output.Time	= 0.;
		leaf_output.Good	= -9999.;
	
		leaf_output.InTime	= 0;
		leaf_output.Wall	= 0;
		
		for ( int iType = 0; iType < 3; iType++ ) {
			leaf_output.n50		[iType]		= 0;
			leaf_output.dir		[iType][0]	= 0;
			leaf_output.dir		[iType][1]	= 0;
			leaf_output.dir		[iType][2]	= 0;
			leaf_output.dir_goodness[iType]		= 0;
			leaf_output.dirKS	[iType]		= 0;
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

		/*bool bID =*/ AnalyseEvent(fIDevent,ID_EVENT);

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
		/*bool bmPMT =*/ AnalyseEvent(fHybridevent,mPMT_EVENT);
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
		/*bool bOD =*/ AnalyseODEvent(fODevent,OD_EVENT);

		Hit_OD     = fHit;
		Hit_OD_50  = fHit_50;
		Hit_OD_200 = fHit_200;
		Hit_OD_400 = fHit_400;
#endif	
		/****************************************************************************************/
		/* Benjamin Fitter									*/
		/****************************************************************************************/

		std::cout << " Start BQ " << std::endl;
		TStopwatch timerBQ;
		timerBQ.Reset();
		timerBQ.Start();
		
		// 
		leaf_output = BQFitter::GetME()->MakeFit();
		
		timerBQ.Stop();
		
		fBQTime = timerBQ.RealTime();
		
		std::cout << " n50 " << leaf_output.n50[0] << " " << leaf_output.n50[1] << " " << leaf_output.n50[2] << std::endl;
		std::cout << " BQ took: " << timerBQ.RealTime() << std::endl;

		/****************************************************************************************/
		/* Bonsai										*/
		/****************************************************************************************/
#ifdef BONSAI
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

	fOutputFile->Write("",TObject::kOverwrite);
	
	// Cleaning
	delete fPrimaryTree;
	delete fOutputFile;

	return 1;

}

void SetCustomBranch(TTree* fPrimaryTree) {		

	fPrimaryTree->Branch("eventId",			&eventId,		"eventId/I");	
	fPrimaryTree->Branch("triggerId",		&triggerId,		"triggerId/I");
	
	fPrimaryTree->Branch("true_particleId",        	&true_particleId);
	fPrimaryTree->Branch("true_vertex_X",     	&true_vertex_X);
	fPrimaryTree->Branch("true_vertex_Y",     	&true_vertex_Y);
	fPrimaryTree->Branch("true_vertex_Z",     	&true_vertex_Z);
	fPrimaryTree->Branch("true_vertex_T",     	&true_vertex_T);
	fPrimaryTree->Branch("true_dir_X",     	&true_dir_X);
	fPrimaryTree->Branch("true_dir_Y",     	&true_dir_Y);
	fPrimaryTree->Branch("true_dir_Z",     	&true_dir_Z);

	fPrimaryTree->Branch("true_energy",     	&true_energy);
	
	fPrimaryTree->Branch("rawhit_num", 		&rawhit_num,		"rawhit_num/I");
	fPrimaryTree->Branch("digithit_num", 		&digithit_num,		"digithit_num/I");
	
	fPrimaryTree->Branch("ID_hits", 		&Hit_ID,		"Hit_ID/I");	
	fPrimaryTree->Branch("ID_hits_50", 		&Hit_ID_50,		"Hit_ID_50/I");	
	fPrimaryTree->Branch("ID_hits_200", 		&Hit_ID_200,		"Hit_ID_200/I");	
	fPrimaryTree->Branch("ID_hits_400", 		&Hit_ID_400,		"Hit_ID_400/I");
	
	fPrimaryTree->Branch("mPMT_hits", 		&Hit_mPMT,		"Hit_mPMT/I");	
	fPrimaryTree->Branch("mPMT_hits_50", 		&Hit_mPMT_50,		"Hit_mPMT_50/I");	
	fPrimaryTree->Branch("mPMT_hits_200", 		&Hit_mPMT_200,		"Hit_mPMT_200/I");	
	fPrimaryTree->Branch("mPMT_hits_400", 		&Hit_mPMT_400,		"Hit_mPMT_400/I");
	
	fPrimaryTree->Branch("bs_vertex", 		bs_vertex,		"bs_vertex[4]/F");	
	fPrimaryTree->Branch("bs_good", 		bs_good,		"bs_good[3]/F");	
	fPrimaryTree->Branch("bs_ctime", 		&fBSTime,		"bs_ctime/F"); // Computation time
	
	fPrimaryTree->Branch("lf_vertex", 		&leaf_output.Vtx,		"lf_vertex[3]/D");
	fPrimaryTree->Branch("lf_time", 		&leaf_output.Time,		"lf_time/D");
	fPrimaryTree->Branch("lf_intime", 		&leaf_output.InTime,		"lf_intime/I");
	fPrimaryTree->Branch("lf_good", 		&leaf_output.Good,		"lf_good/D");
	fPrimaryTree->Branch("lf_wall", 		&leaf_output.Wall,		"lf_wall/D");
	fPrimaryTree->Branch("lf_n50", 			&leaf_output.n50,		"lf_n50[3]/I");
	fPrimaryTree->Branch("lf_dir", 			&leaf_output.dir,		"lf_dir[3][3]/D");
	fPrimaryTree->Branch("lf_dir_goodness", 	&leaf_output.dir_goodness,	"lf_dir_goodness[3]/D");
	fPrimaryTree->Branch("lf_dirKS", 		&leaf_output.dirKS,		"lf_dirKS[3]/D");
	fPrimaryTree->Branch("lf_ctime", 		&fBQTime,			"lf_ctime/F"); // Computation time

}

void SetGeoBranch(TTree* fGeoTree) {	

	fGeoTree->Branch("wcsimrootgeom",		fGeometry);
	
}

bool AnalyseEvent(WCSimRootEvent * tEvent, int iEventType) {

	int iHybrid = 0;
	if ( iEventType == mPMT_EVENT ) iHybrid = 1;
	
	int   nVertex			= 0; // Number of Vertex in event
	int   nTrack			= 0; // Number of Track in event
	int   nRawCherenkovHits		= 0; // Number of Raw Cherenkov hits
	int   nDigitizedCherenkovHits	= 0; // Number of Digitized Cherenkov hits
	double 	fTriggerTime		= 0.;

	// Declare some useful variables
	WCSimRootTrigger * 	fRootTrigger;
	//TClonesArray *		fTimeArray;
	
	// Currently only one Trigger is used (to be check)
	int iTrig = 0;
	//for(int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++){

		triggerId    = iTrig;
		fRootTrigger = tEvent->GetTrigger(iTrig);
				
		// Grab the big arrays of times and parent IDs
		//fTimeArray   = fRootTrigger->GetCherenkovHitTimes();
			
		// Get number of vertex and tracks
#ifdef OLD_WCSIM
		nVertex                 = 1;
#else
		nVertex                 = fRootTrigger->GetNvtxs();
#endif
		nTrack                  = fRootTrigger->GetNtrack();
			
		if ( nTrack == 0 ) {
			// No track, no hit, nothing to do
			return false;
		}
			
		// Get number of hits
		nRawCherenkovHits	= fRootTrigger->GetNumTubesHit();
		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();
		if((fRootTrigger->GetTriggerInfo()).size()>=3) fTriggerTime            = fRootTrigger->GetTriggerInfo()[2] - fRootTrigger->GetTriggerInfo()[1];
		else{
		  std::cout<<"NO TRIGGER. CHEK THIS"<<std::endl;
		  fTriggerTime=0;
		}
			
		// Loop on vertex
		for(int iVertex = 0; iVertex < nVertex; iVertex++){
#ifdef OLD_WCSIM
		  true_vertex_X.push_back(fRootTrigger->GetVtx(0));
		  true_vertex_Y.push_back(fRootTrigger->GetVtx(1));
		  true_vertex_Z.push_back(fRootTrigger->GetVtx(2));
#ifdef VERBOSE
		  std::cout<<"True vertex position = "<<fRootTrigger->GetVtx(0)<<", "<<fRootTrigger->GetVtx(1)<<", "<<fRootTrigger->GetVtx(2)<<std::endl;
#endif
#else
		  true_vertex_X.push_back(fRootTrigger->GetVtxs(iVertex,0));
		  true_vertex_Y.push_back(fRootTrigger->GetVtxs(iVertex,1));
		  true_vertex_Z.push_back(fRootTrigger->GetVtxs(iVertex,2));
#endif
			if ( iVertex*2 > nTrack ) {
				std::cout << "ERROR: Vertex and Track number incompatible (nVertex: " << nVertex << " ; nTrack " << nTrack << ") " << std::endl;
				return 0;
			}
				
			// Beam info are registered in the Track
			TObject * Track = (fRootTrigger->GetTracks())->At(iVertex*2);
			WCSimRootTrack *wcTrack = dynamic_cast<WCSimRootTrack*>(Track);
				
			true_particleId.push_back( wcTrack->GetIpnu() );
			true_energy.push_back(	   wcTrack->GetE()    );
			true_vertex_T.push_back( wcTrack->GetTime() );
#ifdef V2BQDIR
			// Loop through elements in the TClonesArray of WCSimTracks
			for (int t=0; t<nTrack; t++)
			  {
			    TObject *element = (fRootTrigger->GetTracks())->At(t);

			    //WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
			    WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);

			    if(t==(nTrack-1)){//Mother particle
			      true_dir_X.push_back(wcsimroottrack->GetDir(0));
			      true_dir_Y.push_back(wcsimroottrack->GetDir(1));
			      true_dir_Z.push_back(wcsimroottrack->GetDir(2));
			    }
			  }
#endif
			      
		}
						
		// Send True Position to LEAF for check
		std::vector<double> vVtxTrue(3,0.);
		vVtxTrue[0] = true_vertex_X[0];
		vVtxTrue[1] = true_vertex_Y[0];
		vVtxTrue[2] = true_vertex_Z[0];

		BQFitter::GetME()->SetTrueVertexInfo(vVtxTrue,0);
		
		/*
		std::vector<Hit> RawHit;
		// Loop on Raw Hits	
		for(int iHit = 0; iHit < nRawCherenkovHits; iHit++){
			
			TObject *Hit = (fRootTrigger->GetCherenkovHits())->At(iHit);
			
			WCSimRootCherenkovHit *wcHit = dynamic_cast<WCSimRootCherenkovHit*>(Hit);

			int pmtId      = wcHit->GetTubeID();
				
			int HitTimeIdx = wcHit->GetTotalPe(0);
			int HitPE      = wcHit->GetTotalPe(1);
				
			for (int iHitTime = HitTimeIdx; iHitTime < (HitTimeIdx + HitPE); iHitTime++) {
				// Look for all Hits on one PMT
				WCSimRootCherenkovHitTime * wcHitTime = dynamic_cast<WCSimRootCherenkovHitTime*>(fTimeArray->At(iHitTime));
				// How can we use wcHitTime->GetParentID() ?
				
				struct Hit hHit;
					
				hHit.pmtId = pmtId;
				hHit.T     = (float) wcHitTime->GetTruetime();
				hHit.Q     = 1.;
				hHit.dark  = false;
				hHit.Type  = iEventType;
					
				RawHit.push_back(hHit);
			}
				
			fHit += 1;
		} 
		*/
		if ( nRawCherenkovHits < 1 ) return false;
		
		/*		
		// Sort by time
		std::sort(RawHit.begin(), RawHit.end(), SortingT());
		
		// Get first hit for normalization
		float fFirstHitTime = RawHit[0].T;
		fLastRawHit = RawHit.size();
		*/

		std::vector<float> times;
		
		// Get number of hit
		rawhit_num = nRawCherenkovHits;
		/*	
		for (int iHit = 0; iHit < rawhit_num; iHit++ ) {
			
			Hit hHit = RawHit[iHit];	
							
			//hHit.T -= fFirstHitTime;
						
			rawhit_pmtId.push_back(hHit.pmtId);  		
			rawhit_T    .push_back(hHit.T    );  		
			rawhit_Q    .push_back(hHit.Q    );  		
			rawhit_dark .push_back(hHit.dark );  
		}
		*/
			
		// Loop on Digitized Hits
		int nDigitizedCherenkovHits_before = digithit_pmtId.size();
		for(int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
			
			TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
			WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);

			int pmtId      = wcDigitHit->GetTubeId();			
			float HitT     = wcDigitHit->GetT() + fTriggerTime;
			float HitQ     = wcDigitHit->GetQ();
				
			digithit_pmtId.push_back(pmtId);
			digithit_T.push_back(HitT);
			digithit_Q.push_back(HitQ);
			digithit_dark.push_back(false);
			digithit_Type.push_back(iEventType);
		}
			
		int iIdx_BS = 0;
		
		//std::cout<<"Number of digi hits = "<<nDigitizedCherenkovHits<<" vs "<<digithit_pmtId.size()<<std::endl;
		int nDigitizedCherenkovHits_after = digithit_pmtId.size();
		digithit_num = nDigitizedCherenkovHits_after;
			
		// Feed fitter
		for(int iDigitHit = nDigitizedCherenkovHits_before; iDigitHit < nDigitizedCherenkovHits_after; iDigitHit++){//Modified by B.Q on 2020/05/25 to avoid  double counting in multi-PMT option
					
			times.push_back(digithit_T[iDigitHit]);

			//std::cout << " Add Hit with " << digithit_pmtId[iDigitHit] << " Hybrid: " << iHybrid << std::endl;
			if ( digithit_pmtId[iDigitHit] <= 0 ) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
			BQFitter::GetME()->AddHit(digithit_T[iDigitHit],digithit_Q[iDigitHit],iHybrid,digithit_pmtId[iDigitHit]);
			
			// Bonsai (do not store mPMT hits)
			if ( iIdx_BS < 2000 && digithit_T[iDigitHit] < 4000. && digithit_T[iDigitHit] > -4000. && iHybrid == 0 ) {
				//std::cout << pmtId << " " << HitT << std::endl;
				bsCAB[iIdx_BS] = digithit_pmtId[iDigitHit];
				bsT  [iIdx_BS] = digithit_T[iDigitHit]; // shift BS time is needed if interaction time is < 0, this needs to be considered
				bsQ  [iIdx_BS] = digithit_Q[iDigitHit];
					
				iIdx_BS += 1;
			}
		}
				
		// Bonsai
		bsnhit[0]               = iIdx_BS;
		fHit			= digithit_pmtId.size();	
					
		// Compute hits
		std::vector<float> Hit_time_50;
		std::vector<float> Hit_time_200;
		std::vector<float> Hit_time_400;
			
		std::sort(times.begin(),times.end());
		for(int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++){
			
			float HitT = times[iDigitHit];
				
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
	
	//}	


	return true;
}

