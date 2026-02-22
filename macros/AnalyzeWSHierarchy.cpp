#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using namespace std;


bool separatedTriggers=false;   //Assume two independent triggers, one for mPMT, one for B&L
bool produceDirectionPDF=true;  
bool thetaOnly=false;
WCSimRootGeom *geo = 0; 

bool eventDisplay_closeWall = false;
int stopCounterPlots=0;
int nPlotsStopping = 50;
int topPMT=19;
int nPMTpermPMT=19;

// Simple example of reading a generated Root file
const int nPMTtypes = 2;
double PMTradius[nPMTtypes];
const int nGroupsPMTs = 3;  //Different PMT config in an mPMT. I assumed here a rotational symetry of the mPMT, so PMT 1 to 12 are the same, 13 to 18 are the same and 19 is separated

// TOF residuals
TH1D * TimeProfile[nPMTtypes];
TH1D * TimeHitProfile[nPMTtypes];
TH1D * TimeTOFProfile[nPMTtypes];
TH1D * HitTimeTOFProfile[nPMTtypes];

// Cherenkov angle
TH1D * ChargeProfile[nPMTtypes];
TH1D * HitProfile[nPMTtypes];

// Charge
TH1D * TotalCharge[nPMTtypes];
TH1D * TotalHit[nPMTtypes];


/*****************************************************************************************************/
/* DEFINE STRUCT AND USEFUL FUNCTIONS */
/*****************************************************************************************************/


struct arguments
{
  char * inputFile  = NULL;
  char * outputFile = NULL;
  int startEvent    = 0;		// First event to analyze
  int endEvent      = 0;		// Last event to analyze
  bool HK           = false;
  bool verbose      = false;
};

arguments FetchInput(int argc, char* argv[])
{
    int c = -1;
    arguments arglist;

    //Input in c the argument (-f etc...) and in optarg the next argument
    //When the above test becomes -1, it means it fails to find a new argument
    std::cout << std::endl;
    while( (c = getopt(argc, argv, "i:o:s:e:h:v")) != -1 )
    {
        switch(c)
        {
        //Input file name
        case 'i':
            arglist.inputFile = optarg;
            std::cout << "Input WCSim file: " << arglist.inputFile << std::endl;
            break; 

        //Output file name
        case 'o':
            arglist.outputFile = optarg;
            if(arglist.outputFile == NULL){sprintf(arglist.outputFile,"out.txt");}
            std::cout << "Output root file: " << arglist.outputFile << std::endl;
            break;

        //Starting event
        case 's':
            arglist.startEvent = atoi(optarg);
            if(arglist.startEvent == NULL  ||  arglist.startEvent < 0){arglist.startEvent = 0;}
            std::cout << "Starting event #" << arglist.startEvent << std::endl;
            break;

        //Ending event
        case 'e':
            arglist.endEvent = atoi(optarg);
            if(arglist.endEvent == NULL){arglist.endEvent = 0;}
			if(arglist.endEvent >= arglist.startEvent){std::cout << "Ending event #" << arglist.endEvent << std::endl;}
			if(arglist.endEvent <  arglist.startEvent){std::cout << "Ending event = last WCSim event" << std::endl;}
            break;

        //Using HK geometry
        case 'h':
            arglist.HK = true;
            std::cout << "Using HK geoemtry" << std::endl;
            break;

        //Warning
        case 'v':
            arglist.verbose = true;
            std::cout << "VERBOSE option on" << std::endl;
            break;
        }
    }
    std::cout << std::endl;

    return arglist;
}


/*****************************************************************************************************/
/* MAKE PLOTS */
/*****************************************************************************************************/

int main(int argc, char **argv)
{
	//Get arguments
  arguments arglist = FetchInput(argc, argv);
  char * filename    = arglist.inputFile;
  char * outfilename = arglist.outputFile;
  int startEvent = arglist.startEvent;	//We start counting at 0
  int endEvent   = arglist.endEvent;
  bool HK        = arglist.HK;
  bool verbose   = arglist.verbose;

  bool hybrid = true;

  float cvacuum = 0.299792458;//speed of light, in meter per ns.
  float nindex = 1.373;//refraction index of water
  double TotalNHits = 1e4;
  double TankSize = 1100/2;//in cm, the maximal size of the plots for vertices etc...
  double TankRadius = 742/2;//in cm, the maximal size of the plots for vertices etc...
  double TankHalfHeight = 1042/2;//in cm, the maximal size of the plots for vertices etc...
  if(HK)
  {
    TotalNHits = 3e4;
    TankSize = 7100;//in cm, the maximal size of the plots for vertices etc...
    TankRadius = 3242.766;//in cm, the maximal size of the plots for vertices etc...
    TankHalfHeight = 3296.471;//in cm, the maximal size of the plots for vertices etc...
  }




  // Open the input file
  TFile *file;
  if (filename==NULL)
  {
    cout << "Error, no input file" << endl;
    return -1;
  }
  else
  {
    file = new TFile(filename,"read");
  }

  if (!file->IsOpen())
  {
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }

  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  int nevent = ((int)tree->GetEntries());
  if(endEvent>startEvent  &&  endEvent!=0  &&  endEvent<=nevent) nevent = endEvent;
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent  = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);    // Force deletion to prevent memory leak 

  TBranch *branch2;
  if(hybrid)
  {
    branch2 = tree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
    tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);     // Force deletion to prevent memory leak 
  }

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {exit(9);}
  geotree->GetEntry(0);
  PMTradius[0]=geo->GetWCPMTRadius();
  PMTradius[1]=geo->GetWCPMTRadius(true);
  cout << "Number of PMTs of 1st type = " << geo->GetWCNumPMT()     << ", radius = " << PMTradius[0] << endl;
  cout << "Number of PMTs of 2nd type = " << geo->GetWCNumPMT(true) << ", radius = " << PMTradius[1] << endl;
  
  // Options tree - only need 1 "event"
  TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
  WCSimRootOptions *opt = 0; 
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) 
  {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // Start with the main "subevent", as it contains most of the info and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;

  TH1F *hvtx0 = new TH1F("Event VTX0", "Event VTX0", 200, -TankRadius, TankRadius);
  TH1F *hvtx1 = new TH1F("Event VTX1", "Event VTX1", 200, -TankRadius, TankRadius);
  TH1F *hvtx2 = new TH1F("Event VTX2", "Event VTX2", 200, -TankHalfHeight, TankHalfHeight);
  
  int num_trig=0;
  int nbins_position = TankSize/2;
  int nbins_angle = 720;
  int nbins_time = 1e4;
  int nbins_TOF = 1e3;
  int nbins_pe = 500;
  int nbins_totalcharge = 3e4;
  int nbins_distancePMT = 700;
  int nbins_phi = 18;
  int nbins_theta = 18;
  int nbins_weight = 200;  

  for(int i=0 ; i<nPMTtypes ; i++)
  {
    // TOF residuals
    TimeProfile[i]      = new TH1D(Form("TimeProfile_pmtType%d",i),"",      nbins_time, 0, 5e3);          TimeProfile[i]->Sumw2();
    TimeHitProfile[i]   = new TH1D(Form("TimeHitProfile_pmtType%d",i),"",   nbins_time, 0, 5e3);          TimeHitProfile[i]->Sumw2();
    TimeTOFProfile[i]   = new TH1D(Form("TimeTOFProfile_pmtType%d",i),"",   nbins_time, -1.5e3, 1.5e3);   TimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfile[i]= new TH1D(Form("HitTimeTOFProfile_pmtType%d",i),"",nbins_time, -1.5e3, 1.5e3);   HitTimeTOFProfile[i]->Sumw2();

    // Cherenkov angle
    ChargeProfile[i] = new TH1D(Form("ChargeProfile_pmtType%d",i),"",nbins_angle, 0, 180);                ChargeProfile[i]->Sumw2(); 
    HitProfile[i]    = new TH1D(Form("HitProfile_pmtType%d",i),"",   nbins_angle, 0, 180);                HitProfile[i]->Sumw2();

    // Charge
    TotalCharge[i]  = new TH1D(Form("TotalCharge_pmtType%d",i),"",  nbins_totalcharge, 0, 3e4);           TotalCharge[i]->Sumw2();
    TotalHit[i]     = new TH1D(Form("TotalHit_pmtType%d",i),"",     nbins_totalcharge, 0, TotalNHits);    TotalHit[i]->Sumw2();
  }



  /**************/
  /* LOOP EVENT */
  /**************/

  for (int ev=startEvent ; ev<nevent ; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);

    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);

    std::vector<double> triggerInfo;              triggerInfo  = wcsimrootevent->GetTriggerInfo();
    std::vector<double> triggerInfo2;  if(hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();

    // Verbose print
    if(verbose)
    {
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(), wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n", wcsimrootsuperevent->GetNumberOfSubEvents());
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0), wcsimrootevent->GetVtx(1), wcsimrootevent->GetVtx(2));
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());

      for(int v=0 ; v<triggerInfo.size() ; v++)            {cout << "Trigger entry #"  << v << ", info = " << triggerInfo[v]  << endl;}
      if(hybrid) for(int v=0 ; v<triggerInfo2.size() ; v++){cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << endl;}
    }

    // Get Trigger info
    double triggerShift[nPMTtypes];
    double triggerTime[nPMTtypes];
    for(int pmtType=0 ; pmtType<nPMTtypes ; pmtType++)
    {
      triggerShift[pmtType] = 0;
      triggerTime[pmtType]  = 0;

      if(triggerInfo.size() >= 3)
      {
        if(pmtType==0)
        {
          triggerShift[pmtType] = triggerInfo[1];
          triggerTime[pmtType]  = triggerInfo[2];
        }
      }

      if(triggerInfo2.size() >= 3)
      {
        if(pmtType==1 && hybrid)
        {
          triggerShift[pmtType] = triggerInfo2[1];
          triggerTime[pmtType]  = triggerInfo2[2];
        }
      }
    }



    /**************/
    /* LOOP TRACK */
    /**************/

    int ntrack = wcsimrootevent->GetNtrack();    // Get the number of tracks
    double particleStart[3];
    double particleStop[3];
    double particleDir[3];

    // Loop through elements in the TClonesArray of WCSimTracks
    for (int i=0 ; i<ntrack ; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);

      // Mother particle
      if(i == (ntrack-1))
      {
        for(int j=0; j<3; j++)
        {
          particleStart[j] = wcsimroottrack->GetStart(j);
          particleStop[j]  = wcsimroottrack->GetStop(j);
          particleDir[j]   = wcsimroottrack->GetDir(j);
        }
      }
    }



    /*************/
    /* LOOP HITS */
    /*************/

    // Get the number of Cherenkov hits.
    // Note... this is *NOT* the number of photons that hit tubes.
    // It is the number of tubes hit with Cherenkov photons.
    // The number of digitized tubes will be smaller because of the threshold.
    // Each hit "raw" tube has several photon hits.  The times are recorded.
    // For digitized info (one time/charge tube after a trigger) use the digitized information.

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovhits2 = 0;     
    if(hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();

    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
    int ncherenkovdigihits2 = 0; 
    if(hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits(); 

    // Verbose print
    if(verbose)
    {
      cout << "LOOP OVER DIGITIZED HITS:" << endl;
      printf("node id: %i\n", ev);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
    }

    // Loop over PMT type
    for(int pmtType=0 ; pmtType<nPMTtypes ; pmtType++)
    {
      if(separatedTriggers)
      {
        if(triggerInfo2.size()!=0 && pmtType==0) continue;
        if(triggerInfo.size()!=0  && pmtType==1) continue;
      }

      int nhits;
      if(pmtType == 0) nhits = ncherenkovdigihits;
      else             nhits = ncherenkovdigihits2;

      for (int i=0 ; i<nhits ; i++)
      {
        // Get PMT info
        TObject *Hit;
        if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
        else           Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

        WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
        int tubeNumber           = wcsimrootcherenkovdigihit->GetTubeId();
        double peForTube         = wcsimrootcherenkovdigihit->GetQ();

        WCSimRootPMT pmt;
        if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1, false);
        else             pmt = geo->GetPMT(tubeNumber-1, true); 
        
        double PMTpos[3];
        for(int j=0 ; j<3 ; j++){PMTpos[j] = pmt.GetPosition(j);}


        // Compute variables
        double vDir[3];
        for(int j=0 ; j<3 ; j++){vDir[j] = PMTpos[j] - particleStart[j];}

        double Norm = TMath::Sqrt(vDir[0]*vDir[0] + vDir[1]*vDir[1] + vDir[2]*vDir[2]);
        for(int j=0 ; j<3 ; j++){vDir[j]   /= Norm;}

        double time = wcsimrootcherenkovdigihit->GetT();
        double tof = Norm*1e-2/(cvacuum/nindex);
        double tof_residuals = time-tof+triggerTime[pmtType]-triggerShift[pmtType];
        double relativeAngle = TMath::ACos(vDir[0]*particleDir[0] + vDir[1]*particleDir[1] + vDir[2]*particleDir[2])*180./TMath::Pi();


        // Fill histograms
        // TOF residuals
        TimeProfile[pmtType]->Fill(time, peForTube);
        TimeHitProfile[pmtType]->Fill(time, 1.);
        TimeTOFProfile[pmtType]->Fill(tof_residuals, peForTube);
        HitTimeTOFProfile[pmtType]->Fill(tof_residuals, 1.);

        // Cherenkov angle
        ChargeProfile[pmtType]->Fill(relativeAngle, peForTube);
        HitProfile[pmtType]->Fill(relativeAngle, 1);

        // Charge
        TotalCharge[pmtType]->Fill(peForTube);
        TotalHit[pmtType]->Fill(1.);
      }
    }
    // Vertex
    hvtx0->Fill(wcsimrootevent->GetVtx(0));
    hvtx1->Fill(wcsimrootevent->GetVtx(1));
    hvtx2->Fill(wcsimrootevent->GetVtx(2));



    /****************/
    /* REINITIALIZE */
    /****************/

    wcsimrootsuperevent->ReInitialize();
    if(hybrid) wcsimrootsuperevent2->ReInitialize();
  }



  // Write to file
  cout << "num_trig " << num_trig << "\n";  
  TFile * outfile = new TFile(outfilename,"RECREATE");
  cout << "File " << outfilename << " is open for writing" << endl;
  for(int i=0 ; i<nPMTtypes ; i++)
  {
    // TOF residuals
    TimeProfile[i]->Write();
    TimeHitProfile[i]->Write();
    TimeTOFProfile[i]->Write();
    HitTimeTOFProfile[i]->Write();

    // Cherenkov angle
    ChargeProfile[i]->Write();
    HitProfile[i]->Write();

    // Charge
    TotalCharge[i]->Write();
    TotalHit[i]->Write();
  }

  // Vertex
  hvtx0->Write();
  hvtx1->Write();
  hvtx2->Write();
  
  outfile->Close();
  
  return 0;
}