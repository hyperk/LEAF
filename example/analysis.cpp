/*********************************************************************************/
/**	analysis.cpp								**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		**/
/**	Date: Febuary 10th 2020							**/
/**	Desc: Example application code for Benjamin's Low-E Fitter for Hyper-K	**/
/*********************************************************************************/

#include "analysis.hpp"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///* Some Utility Functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double GetResidualTime(const std::vector<double>& origin, double originTime, WCSimRootCherenkovDigiHit *wcDigitHit, TimeDelta fTimeCorrection, double fLfTriggerTime)
{
	float fCVacuum = 29.9792458; // speed of light, in centimeter per ns.
	float fNIndex = 1.373;			  // 1.385;//1.373;//refraction index of water
	double fLightSpeed = fCVacuum / fNIndex;
	double HitT = wcDigitHit->GetT() + fLfTriggerTime;
	WCSimRootPMT pmt = fLeafGeometry->GetPMT(wcDigitHit->GetTubeId() - 1, false);
	double PMTpos[3];
	for (int j = 0; j < 3; j++) PMTpos[j] = pmt.GetPosition(j);
	double distance = Astro_GetDistance(PMTpos, origin);
	double tof = distance / fLightSpeed;
	double hitTime = (fTimeCorrection + HitT) / TimeDelta::ns;
	return hitTime - tof - originTime;
}

arguments FetchInput(int argc, char** argv)
{
    int c = -1;
    arguments arglist;

    //Input in c the argument (-f etc...) and in optarg the next argument
    //When the above test becomes -1, it means it fails to find a new argument
    std::cout << std::endl;
    while( (c = getopt(argc, argv, "i:o:d:h:s:e:v")) != -1 )
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
            if(arglist.outputFile == ""){ arglist.outputFile = "out.txt";}
            std::cout << "Output root file: " << arglist.outputFile << std::endl;
            break;

        //Darknoise
        case 'd':
            arglist.darkNoise = atof(optarg);
            std::cout << "Dark noise frequency: " << arglist.darkNoise << " kHz" << std::endl;
            break;

        //Darknoise hybrid
        case 'h':
            arglist.darkNoiseH = atof(optarg);
            arglist.hybrid     = true;
            std::cout << "Dark noise frequency (hybrid geometry): " << arglist.darkNoiseH << " kHz" << std::endl;
            break;

        //Starting event
        case 's':
            arglist.startEvent = atoi(optarg);
            if(arglist.startEvent < 0){arglist.startEvent = 0;}
            std::cout << "Starting event #" << arglist.startEvent << std::endl;
            break;

        //Ending event
        case 'e':
            arglist.endEvent = atoi(optarg);
            if(arglist.endEvent < 0){arglist.endEvent = 0;}
			if(arglist.endEvent >= arglist.startEvent){std::cout << "Ending event #" << arglist.endEvent << std::endl;}
			if(arglist.endEvent <  arglist.startEvent){std::cout << "Ending event = last WCSim event" << std::endl;}
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* MAIN EXECUTION */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	//Get arguments
    arguments arglist = FetchInput(argc, argv);
    std::string sInputFile  = arglist.inputFile;
    std::string sOutputFile = arglist.outputFile;
    double dDarkNoise       = arglist.darkNoise;
    double dDarkNoiseHybrid = 0.;
#ifdef mPMT
	dDarkNoiseHybrid = arglist.darkNoiseH;
#endif
	int firstEvents = arglist.startEvent;	//We start counting at 0
	int lastEvents  = arglist.endEvent;
    bool verbose    = arglist.verbose;

	// Read WCSim output
	TFile *fInputFile = new TFile(sInputFile.c_str(), "READ");

	// Get TTrees
	TTree *fInputTree    = (TTree *)fInputFile->Get("wcsimT");
	TTree *fInputGeoTree = (TTree *)fInputFile->Get("wcsimGeoT");

	fLeafGeometry = 0;
	fInputGeoTree->SetBranchAddress("wcsimrootgeom", &fLeafGeometry);

	// Set Branch for PMTs
	WCSimRootEvent *fIDevent = new WCSimRootEvent();
	fInputTree->SetBranchAddress("wcsimrootevent", &fIDevent);
	fInputTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);	    // Set autodelete to avoid memory leak

	// Set Branch for mPMTs
#ifdef mPMT
	WCSimRootEvent *fHybridevent = new WCSimRootEvent();
	fInputTree->SetBranchAddress("wcsimrootevent2", &fHybridevent);
	fInputTree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);	    // Set autodelete to avoid memory leak
#endif

	// Set Branch for OD PMTs
#ifdef OD_ON
	WCSimRootEvent *fODevent = new WCSimRootEvent();
	fInputTree->SetBranchAddress("wcsimrootevent_OD", &fODevent);
	fInputTree->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);	// Set autodelete to avoid memory leak
#endif

	// Create Output TTree
	TFile *fOutputFile = new TFile(sOutputFile.c_str(), "RECREATE");
	fOutputFile->SetCompressionLevel(2);

	// Read Geo and set branch
	fInputGeoTree->GetEntry(0);
	TTree *fGeoTree     = new TTree("wcsimGeoT", "Geometry TTree");
	TTree *fPrimaryTree = new TTree("Reduced", "Reduced TTree");
	fGeoTree->Branch("wcsimrootgeom", fLeafGeometry);
	SetCustomBranch(fPrimaryTree);
	fGeoTree->Fill();

	// Get PMT Number:
	int nPMT_ID = fLeafGeometry->GetWCNumPMT();
	std::cout << " ID " << nPMT_ID << std::endl;
#ifdef OD_ON
	int nPMT_OD = fLeafGeometry->GetODWCNumPMT();
	std::cout << " OD " << nPMT_OD << std::endl;
#endif
	int nMultPMT = fLeafGeometry->GetWCNumPMT(true);
	std::cout << " mPMT " << nMultPMT << std::endl;

	// Initialize HK Manager geometry
	HKManager::GetME()->SetGeometry(fLeafGeometry, dDarkNoise * 1e3, dDarkNoiseHybrid * 1e3);

	// Initialize LEAF
	LEAF::GetME()->Initialize(HKManager::GetME()->GetGeometry()); // Loads the geometry in leaf and loads the pdfs

	//* for each config variable, check if it has been set, take the default variable if not
	const char* envValue = std::getenv("maxHitsAngle");
	if (envValue != nullptr) 
	{
		maxHitAngle = std::stof(envValue);
		std::cout << "max hit angle : " << maxHitAngle << "°" << std::endl;
	} else std::cout << "Environment variable max hit angle is not set. Using default value : " << maxHitAngle << std::endl;

	envValue = std::getenv("N_Neighbors");
	if (envValue != nullptr) 
	{
		N_Neighbors = std::stof(envValue);
		std::cout << "N_Neighbors : " << N_Neighbors << std::endl;
	} else std::cout << "Environment variable N_Neighbors is not set. Using default value : " << N_Neighbors << std::endl;

	envValue = std::getenv("maxDistanceToNeighbors");
	if (envValue != nullptr) 
	{
		maxDistanceToNeighbors = std::stof(envValue);
		std::cout << "maxDistanceToNeighbors : " << maxDistanceToNeighbors << "cm" << std::endl;
	} else std::cout << "Environment variable maxDistanceToNeighbors is not set. Using default value : " << maxDistanceToNeighbors << std::endl;

	// Read Input Tree
	int iWrite = 0;
	int failAmount = 0;
	int nPrimaryEvents = fInputTree->GetEntries();
	if(lastEvents > 0) nPrimaryEvents = std::min(nPrimaryEvents, lastEvents);

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	// Loop on Primary events
	for (int i=firstEvents ; i<nPrimaryEvents ; i++)
	{
		// Reset Hit vector and clock
		HKManager::GetME()->ResetHitInfo();
		// HKManager::GetME()->ResetSecondaryHitInfo();
		timer.Stop();
		timer.Reset();

		if(verbose)
		{
			std::cout << "\n===========================================================================================================================================================" << std::endl;
			std::cout << "Event # = " << i+1 << " / " << nPrimaryEvents << std::endl;
		}

		timer.Start();
		fInputTree->GetEntry(i);

		// Initialize output variables
		eventId = iWrite;
		nTrigger = 0;
		usedTriggerId = 0;
		rawhit_num_noDN = 0;

		true_particleId.clear();
		true_energy.clear();
		true_origin_X.clear();
		true_origin_Y.clear();
		true_origin_Z.clear();
		true_origin_T.clear();

		leaf_Vertex.clear();
		leaf_Dir.clear();
		leaf_MyDir.clear();
		leaf_QuickDir.clear();

		digithit_pmtId.clear();
		digithit_T.clear();
		correctedDigithit_T.clear();
		digithit_Q.clear();
		Charge_PMT.clear();

		rawhit_num = 0;
		digithit_num = 0;

		fLFTime = -9999;

		Hit_ID = 0;
		Hit_ID_50 = 0;
		Hit_ID_200 = 0;
		Hit_ID_400 = 0;

		Hit_mPMT = 0;
		Hit_mPMT_20 = 0;
		Hit_mPMT_50 = 0;
		Hit_mPMT_200 = 0;
		Hit_mPMT_400 = 0;

		Hit_OD = 0;
		Hit_OD_50 = 0;
		Hit_OD_200 = 0;
		Hit_OD_400 = 0;

		fLastRawHit = 0;

		trueVertex = std::vector<double>(4);
		trueDir = std::vector<double>(3);

		bestTrigger = 0;
		fLfTriggerTime = 0.;

		/****************************************************************************************/
		/* ID events										*/
		/****************************************************************************************/

		fHit     = 0;
		fHit_20  = 0;
		fHit_50  = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		// this is where we send all the informations of the event to leaf
		// we return false if we don't want to process this event based on defined filters
		bool shouldProcess = AnalyseEvent(fIDevent, ID_EVENT);

		if(!shouldProcess)
		{
			failAmount++;
			std::cout << "Event # = " << i << " cannot be processed " << std::endl;
			continue;
		}

		Hit_ID = fHit;
		Hit_ID_20 = fHit_20;
		Hit_ID_50 = fHit_50;
		Hit_ID_200 = fHit_200;
		Hit_ID_400 = fHit_400;

		/****************************************************************************************/
		/* mPMT events										*/
		/****************************************************************************************/

#ifdef mPMT
		fHit = 0;
		fHit_20 = 0;
		fHit_50 = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		// /*bool bmPMT =*/AnalyseEvent(fHybridevent, mPMT_EVENT);

		Hit_mPMT = fHit;
		Hit_mPMT_20 = fHit_20;
		Hit_mPMT_50 = fHit_50;
		Hit_mPMT_200 = fHit_200;
		Hit_mPMT_400 = fHit_400;
#endif

		/****************************************************************************************/
		/* OD events										*/
		/****************************************************************************************/

#ifdef OD_ON
		fHit = 0;
		fHit_20 = 0;
		fHit_50 = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		// There should be a dedicated OD analyser, as many thing should be different than for ID
		// Doesn't exist yet
		/*bool bOD =*/ AnalyseODEvent(fODevent, OD_EVENT);

		Hit_OD = fHit;
		Hit_OD_50 = fHit_50;
		Hit_OD_200 = fHit_200;
		Hit_OD_400 = fHit_400;
#endif

		/****************************************************************************************/
		/* Benjamin Quilain Fitter									*/
		/****************************************************************************************/

		TStopwatch timerLF;
		timerLF.Reset();
		timerLF.Start();

		// To be replaced by meaningful TriggerTime
		TimeDelta fDummyTrigger(0.);

		fOutputFitter = LEAF::GetME()->MakeSequentialFit(HKManager::GetME()->GetHitCollection(), fDummyTrigger);
		// the fit method can be replaced by a joint fit (doesn't work better for now)
		// fOutputFitter = LEAF::GetME()->MakeJointFit(HKManager::GetME()->GetHitCollection(), fDummyTrigger);
		fOutputProps = LEAF::fOutputProps;
		
		timerLF.Stop();

		fLFTime = timerLF.RealTime();

		if (verbose)
		{
			std::cout << "  LEAF took: " << timerLF.RealTime() << " sec for " << HKManager::GetME()->GetHitCollection()->Size() << " Hits";
			std::cout << "  (Vtx Search : "   << fOutputProps.Vtx_Search_ComputeTime 
			          << " , Vtx Minimize : " << fOutputProps.Vtx_Minimize_ComputeTime 
					  << ", Dir Search : "    << fOutputProps.Dir_Search_ComputeTime 
					  << " , Dir Minimize : " << fOutputProps.Dir_Minimize_ComputeTime 
					  << " , Energy Fit : "   << fOutputProps.Energy_Fit_ComputeTime << ")" << std::endl;
			std::cout << "  True vertex: (" << trueVertex[0]
			                        << ", " << trueVertex[1] 
									<< ", " << trueVertex[2] 
									<< ", " << trueVertex[3] << ")  [cm/cm/cm/ns] "
    			    << "  //  direction: (" << trueDir[0] 
						            << ", " << trueDir[1]
						            << ", " << trueDir[2] 
					    << "  //  energy: " << trueEnergy << std::endl;
			std::cout << "  LEAF vertex: (" << fOutputFitter.Vtx[0] 
			                        << ", " << fOutputFitter.Vtx[1] 
									<< ", " << fOutputFitter.Vtx[2] 
									<< ", " << fOutputFitter.Vtx[3] << ")  [cm/cm/cm/ns]" 
    			    << "  //  direction: (" << fOutputFitter.Dir[0] 
						            << ", " << fOutputFitter.Dir[1]
						            << ", " << fOutputFitter.Dir[2] 
					    << "  //  energy: " << fOutputFitter.Energy << std::endl;
		}

		/****************************************************************************************/
		/* Fill output tree									*/
		/****************************************************************************************/
		
		fPrimaryTree->Fill();
		iWrite += 1;
	}

	std::cout << "\n===========================================================================================================================================================";
	std::cout << "\n===========================================================================================================================================================\n" << std::endl;
	std::cout << ((nPrimaryEvents - failAmount) / nPrimaryEvents) * 100 << "% of the events have been processed" << std::endl;

	fOutputFile->cd(); 
	fOutputFile->Write("", TObject::kOverwrite);

	delete fPrimaryTree;
	delete fOutputFile;
	return 1;
}


//* All the variables that will be kept in the output Tree
void SetCustomBranch(TTree *fPrimaryTree)
{
	fPrimaryTree->Branch("eventId", &eventId, "eventId/I");
	fPrimaryTree->Branch("nTrigger", &nTrigger, "nTrigger/I");
	fPrimaryTree->Branch("triggerUsed", &usedTriggerId, "triggerUsed/I");
	fPrimaryTree->Branch("rawhit_num", &rawhit_num, "rawhit_num/I");
	fPrimaryTree->Branch("digithit_num", &digithit_num, "digithit_num/I");

	fPrimaryTree->Branch("true_particleId", &true_particleId);
	fPrimaryTree->Branch("true_origin_X", &true_origin_X);
	fPrimaryTree->Branch("true_origin_Y", &true_origin_Y);
	fPrimaryTree->Branch("true_origin_Z", &true_origin_Z);
	fPrimaryTree->Branch("true_origin_T", &true_origin_T);
	fPrimaryTree->Branch("true_energy", &true_energy);
	fPrimaryTree->Branch("true_dir", &trueDir, "true_dir[3]/D");

	fPrimaryTree->Branch("ID_hits", &Hit_ID, "Hit_ID/I");
	fPrimaryTree->Branch("ID_hits_50", &Hit_ID_50, "Hit_ID_50/I");
	fPrimaryTree->Branch("ID_hits_200", &Hit_ID_200, "Hit_ID_200/I");
	fPrimaryTree->Branch("ID_hits_400", &Hit_ID_400, "Hit_ID_400/I");

	fPrimaryTree->Branch("mPMT_hits", &Hit_mPMT, "Hit_mPMT/I");
	fPrimaryTree->Branch("mPMT_hits_20", &Hit_mPMT_20, "Hit_mPMT_20/I");
	fPrimaryTree->Branch("mPMT_hits_50", &Hit_mPMT_50, "Hit_mPMT_50/I");
	fPrimaryTree->Branch("mPMT_hits_200", &Hit_mPMT_200, "Hit_mPMT_200/I");
	fPrimaryTree->Branch("mPMT_hits_400", &Hit_mPMT_400, "Hit_mPMT_400/I");

	fPrimaryTree->Branch("lf_vertex", &leaf_Vertex);
	fPrimaryTree->Branch("lf_NLL", &fOutputFitter.NLL, "lf_NLL/D");
	fPrimaryTree->Branch("lf_good", &fOutputFitter.NLLR, "lf_good/D");
	fPrimaryTree->Branch("lf_ctime", &fLFTime, "lf_ctime/D"); // Computation time
	fPrimaryTree->Branch("lf_energy", &fOutputFitter.Energy, "lf_energy/D");
	fPrimaryTree->Branch("lf_Dir", &leaf_Dir);
	fPrimaryTree->Branch("lf_MyDir", &leaf_MyDir);
	fPrimaryTree->Branch("lf_Quick_Dir", &leaf_QuickDir);
	fPrimaryTree->Branch("lf_Dir_NLL", &fOutputFitter.DNLL, "lf_Dir_NLL/D");
	fPrimaryTree->Branch("lf_ComputeTime", &fOutputProps.Leaf_ComputeTime, "Leaf_ComputeTime/D");
	fPrimaryTree->Branch("lf_Vtx_Search_ComputeTime", &fOutputProps.Vtx_Search_ComputeTime, "Vtx_Search_ComputeTime/D");
	fPrimaryTree->Branch("lf_Vtx_Minimize_ComputeTime", &fOutputProps.Vtx_Minimize_ComputeTime, "Vtx_Minimize_ComputeTime/D");
	fPrimaryTree->Branch("lf_Dir_Quick_Search_ComputeTime", &fOutputProps.Dir_Quick_Search_ComputeTime, "Dir_Quick_Search_ComputeTime/D");
	fPrimaryTree->Branch("lf_Dir_Search_ComputeTime", &fOutputProps.Dir_Search_ComputeTime, "Dir_Search_ComputeTime/D");
	fPrimaryTree->Branch("lf_Dir_Minimize_ComputeTime", &fOutputProps.Dir_Minimize_ComputeTime, "Dir_Minimize_ComputeTime/D");
	fPrimaryTree->Branch("lf_Energy_Fit_ComputeTime", &fOutputProps.Energy_Fit_ComputeTime, "Energy_Fit_ComputeTime/D");

	fPrimaryTree->Branch("rawTriggerTime", &rawTriggerTime, "rawTriggerTime/D");
	fPrimaryTree->Branch("TotalCharge", &fOutputFitter.TotalCharge, "TotalCharge/D");
	
	fPrimaryTree->Branch("DigiHitT", &digithit_T);
	fPrimaryTree->Branch("CorrectedDigiHitT", &correctedDigithit_T);
	fPrimaryTree->Branch("HitIsDR", &hit_is_DR);
	fPrimaryTree->Branch("DigiHitQ", &digithit_Q);
	fPrimaryTree->Branch("DigiHitResidual", &hit_residual);
	fPrimaryTree->Branch("Charge_PMT", &Charge_PMT);
}


bool AnalyseEvent(WCSimRootEvent *tEvent, int iEventType)
{
	int iHybrid = 0;
	int nVertex = 0;				 // Number of Vertex in event
	int nTrack = 0;					 // Number of Track in event
	int nRawCherenkovHits = 0;		 // Number of Raw Cherenkov hits
	int nDigitizedCherenkovHits = 0; // Number of Digitized Cherenkov hits
	nTrigger = tEvent->GetNumberOfEvents();
	
	usedTriggerId = -1;
	rawTriggerTime = 0.;
	WCSimRootTrigger *fRootTrigger;

	TimeDelta fDummyTrigger(0.);
	TimeDelta fTimeCorrection = HKManager::GetME()->GetHitCollection()->timestamp - fDummyTrigger; //*still have no idea what this is

	//* find best trigger (best number of hits)
	int bestNumbOfHits = 0;
	bestTrigger = -1;
	for (int iTrig = 0; iTrig < nTrigger; iTrig++)
	{
		fRootTrigger = tEvent->GetTrigger(iTrig);
		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();

		if(nDigitizedCherenkovHits > bestNumbOfHits && fRootTrigger->GetTriggerInfo().size() >= 3)
		{
			bestNumbOfHits = nDigitizedCherenkovHits;
			bestTrigger = iTrig;
		}
	}

	if(bestTrigger < 0) return false; //* Cannot process this event, none of the triggers have sufficent information

	usedTriggerId = bestTrigger;

	//* Retreive true infos (onlyh set by wcsim in the first trigger)

	for (int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++)
	{
		fRootTrigger = tEvent->GetTrigger(iTrig);

#ifdef OLD_WCSIM
		nVertex = 1;
#else
		nVertex = fRootTrigger->GetNvtxs();
#endif

		nTrack = fRootTrigger->GetNtrack();
		nRawCherenkovHits = fRootTrigger->GetNumTubesHit();

		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();

		//* Fist trigger is the one used to compute true Vertex & Dir
		if (iTrig == 0)
		{
			// Loop on vertex
			for (int iVertex = 0; iVertex < nVertex; iVertex++)
			{

#ifdef OLD_WCSIM
				for(int j=0; j<3; j++) trueVertex[j]=fRootTrigger->GetVtx(j);
				true_origin_X.push_back(fRootTrigger->GetVtx(0));
				true_origin_Y.push_back(fRootTrigger->GetVtx(1));
				true_origin_Z.push_back(fRootTrigger->GetVtx(2));
#else
				for(int j=0; j<3; j++) trueVertex[j]=fRootTrigger->GetVtxs(iVertex,j);
				true_origin_X.push_back(fRootTrigger->GetVtxs(iVertex, 0));
				true_origin_Y.push_back(fRootTrigger->GetVtxs(iVertex, 1));
				true_origin_Z.push_back(fRootTrigger->GetVtxs(iVertex, 2));
#endif

				if (iVertex * 2 > nTrack)
				{
					std::cout << "ERROR: Vertex and Track number incompatible (nVertex: " << nVertex << " ; nTrack " << nTrack << ") " << std::endl;
					continue;
				}

				// Beam info are registered in the Track
				TObject *Track = (fRootTrigger->GetTracks())->At(iVertex * 2);
				WCSimRootTrack *wcTrack = dynamic_cast<WCSimRootTrack *>(Track);

				if (wcTrack->GetParentId() == 0 && wcTrack->GetIpnu()==11) 
				{
					for(int j=0; j<3; j++) trueDir[j] = wcTrack->GetPdir(j);  
					Normalize(trueDir);
				}
				
				trueEnergy = wcTrack->GetE();
				true_particleId.push_back(wcTrack->GetIpnu());
				true_energy.push_back(wcTrack->GetE());
				true_origin_T.push_back(fRootTrigger->GetVtx(3));
			}

			// Send True Position to LEAF to make checks of success in the steps of leaf (should be removed for true data)
			std::vector<double> vVtxTrue(3, 0.);
			vVtxTrue[0] = true_origin_X[0];
			vVtxTrue[1] = true_origin_Y[0];
			vVtxTrue[2] = true_origin_Z[0];
			SetTrueVertexInfo(vVtxTrue, true_origin_T[0]);
			SetTrueDirInfo(trueDir);

			rawhit_num = nRawCherenkovHits;
		}

		if(iTrig != usedTriggerId) continue; //* The trigger being processed

		auto triggerInfo = fRootTrigger->GetTriggerInfo();
		fLfTriggerTime = triggerInfo[2] - triggerInfo[1]; // time of th event - offset (common offset of all events, depends on wcsim config giles)
		rawTriggerTime = triggerInfo[2];

		startingCherenkovHitID = digithit_pmtId.size();
		std::vector<double> times;

		// Loop on Digitized Hits
		double AllQ = 0;
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++)
		{
			//* Get Hit Informations
			TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
			WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit *>(Hit);
			
			int pmtId = wcDigitHit->GetTubeId();
			WCSimRootPMT pmt = fLeafGeometry->GetPMT(pmtId - 1, false);
			double HitT   = wcDigitHit->GetT() + fLfTriggerTime;
			double HitQ   = wcDigitHit->GetQ();
			int peForTube = wcDigitHit->GetQ();
			
			//* Fill output root file
			AllQ += HitQ;
			digithit_pmtId.push_back(pmtId);
			digithit_T.push_back(HitT);
			correctedDigithit_T.push_back(HitT - fLfTriggerTime - triggerInfo[1]);
			digithit_Q.push_back(HitQ);
			Charge_PMT.push_back(peForTube);
			digithit_Type.push_back(iEventType);
		}


		int iIdx_BS = 0;

		nDigitizedCherenkovHits = digithit_pmtId.size();
		digithit_num = nDigitizedCherenkovHits - startingCherenkovHitID;

		// Feed fitter
		for (int iDigitHit = startingCherenkovHitID; iDigitHit < nDigitizedCherenkovHits; iDigitHit++)
		{
			times.push_back(digithit_T[iDigitHit]);

			if (digithit_pmtId[iDigitHit] <= 0) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
			HKManager::GetME()->AddHit(digithit_T[iDigitHit], digithit_Q[iDigitHit], iHybrid, digithit_pmtId[iDigitHit]);

			// Bonsai (do not store mPMT hits)
			if (iIdx_BS < 2000 && digithit_T[iDigitHit] < 4000. && digithit_T[iDigitHit] > -4000. && iHybrid == 0)
			{
				// std::cout << pmtId << " " << HitT << std::endl;
				bsCAB[iIdx_BS] = digithit_pmtId[iDigitHit];
				bsT[iIdx_BS] = digithit_T[iDigitHit]; // shift BS time is needed if interaction time is < 0, this needs to be considered
				bsQ[iIdx_BS] = digithit_Q[iDigitHit];

				iIdx_BS += 1;
			}
		}

		// Bonsai
		bsnhit[0] = iIdx_BS;
		fHit = digithit_pmtId.size();

		// Compute hits
		std::vector<double> Hit_time_20;
		std::vector<double> Hit_time_50;
		std::vector<double> Hit_time_200;
		std::vector<double> Hit_time_400;

		std::sort(times.begin(), times.end());
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++)
		{
			double HitT = times[iDigitHit];

			Hit_time_20.push_back(HitT);
			Hit_time_50.push_back(HitT);
			Hit_time_200.push_back(HitT);
			Hit_time_400.push_back(HitT);

			// Count hit in 20 ns window
			while (HitT - Hit_time_20[0] > 20.) Hit_time_20.erase(Hit_time_20.begin());

			// Count hit in 50 ns window
			while (HitT - Hit_time_50[0] > 50.) Hit_time_50.erase(Hit_time_50.begin());

			// Count hit in 200 ns window
			while (HitT - Hit_time_200[0] > 200.) Hit_time_200.erase(Hit_time_200.begin());

			// Count hit in 400 ns window
			while (HitT - Hit_time_400[0] > 400.) Hit_time_400.erase(Hit_time_400.begin());

			if ((unsigned int)fHit_20 < Hit_time_20.size()) fHit_20 = Hit_time_20.size();
			if ((unsigned int)fHit_50 < Hit_time_50.size()) fHit_50 = Hit_time_50.size();
			if ((unsigned int)fHit_200 < Hit_time_200.size()) fHit_200 = Hit_time_200.size();
			if ((unsigned int)fHit_400 < Hit_time_400.size()) fHit_400 = Hit_time_400.size();
		}
	}
	return true;
}
