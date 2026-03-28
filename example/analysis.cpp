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

arguments FetchInput(int argc, char** argv) {
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

int main(int argc, char **argv) {
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
	LEAF::GetME()->Initialize(	HKManager::GetME()->GetGeometry(),
								HKManager::GetME()->GetDarkNoise(),
								HKManager::GetME()->GetGeometryPMT_ID(),
								HKManager::GetME()->GetGeometryPMT_mPMT()); // Loads the geometry in leaf and loads the pdfs

	//* for each config variable, check if it has been set, take the default variable if not
	const char* envValue = std::getenv("maxHitsAngle");
	if (envValue != nullptr) {
		maxHitAngle = std::stof(envValue);
		std::cout << "max hit angle : " << maxHitAngle << "°" << std::endl;
	} 
	else std::cout << "Environment variable max hit angle is not set. Using default value : " << maxHitAngle << std::endl;

	envValue = std::getenv("N_Neighbors");
	if (envValue != nullptr) {
		N_Neighbors = std::stof(envValue);
		std::cout << "N_Neighbors : " << N_Neighbors << std::endl;
	} 
	else std::cout << "Environment variable N_Neighbors is not set. Using default value : " << N_Neighbors << std::endl;

	envValue = std::getenv("maxDistanceToNeighbors");
	if (envValue != nullptr) {
		maxDistanceToNeighbors = std::stof(envValue);
		std::cout << "maxDistanceToNeighbors : " << maxDistanceToNeighbors << "cm" << std::endl;
	} 
	else std::cout << "Environment variable maxDistanceToNeighbors is not set. Using default value : " << maxDistanceToNeighbors << std::endl;

	// Read Input Tree
	int iWrite = 0;
	int failAmount = 0;
	int nPrimaryEvents = fInputTree->GetEntries();
	if(lastEvents > 0) nPrimaryEvents = std::min(nPrimaryEvents, lastEvents);

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	// Loop on Primary events
	for (int i=firstEvents ; i<nPrimaryEvents ; i++) {
		//if ( i > 10 ) break;
		// Reset Hit vector and clock
		HKManager::GetME()->ResetHitInfo();
		// HKManager::GetME()->ResetSecondaryHitInfo();
		timer.Stop();
		timer.Reset();

		if(verbose) {
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
		true_energies.clear();
		true_vertexes.clear();

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
		bool shouldProcess = AnalyseEvent(fIDevent, PMTType::kID);

		if(!shouldProcess) {
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
		/* mPMT events																			*/
		/****************************************************************************************/

#ifdef mPMT
		fHit = 0;
		fHit_20 = 0;
		fHit_50 = 0;
		fHit_200 = 0;
		fHit_400 = 0;

		// /*bool bmPMT =*/AnalyseEvent(fHybridevent, PMTType::kmPMT);

		Hit_mPMT = fHit;
		Hit_mPMT_20 = fHit_20;
		Hit_mPMT_50 = fHit_50;
		Hit_mPMT_200 = fHit_200;
		Hit_mPMT_400 = fHit_400;
#endif

		/****************************************************************************************/
		/* OD events																			*/
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
		/* Benjamin Quilain Fitter																*/
		/****************************************************************************************/

		TStopwatch timerLF;
		timerLF.Reset();
		timerLF.Start();

		LEAF::GetME()->LoadHitsCollection(HKManager::GetME()->GetHitCollection_ID(),HKManager::GetME()->GetHitCollection_mPMT());
		fLeafOutput = LEAF::GetME()->MakeSequentialFit();
		//std::cout << "(X, Y, Z, T): (" << fLeafOutput.vtx.X() << ", " << fLeafOutput.vtx.Y() << ", " << fLeafOutput.vtx.Z() << ", " << fLeafOutput.vtx.T() << ") NLL: " << fLeafOutput.vtx_nll << std::endl;
		//std::cout << "(dX, dY, dZ, T): (" << fLeafOutput.dir.X() << ", " << fLeafOutput.dir.Y() << ", " << fLeafOutput.dir.Z() << ") NLL: " << fLeafOutput.dir_nll << std::endl;

		// the fit method can be replaced by a joint fit (doesn't work better for now)
		// fLeafOutput = LEAF::GetME()->MakeJointFit(HKManager::GetME()->GetHitCollection(), fDummyTrigger);
		fLeafPerformances = LEAF::GetME()->GetPerformances();

		timerLF.Stop();

		fLFTime = timerLF.RealTime();

		if (verbose) {
			std::cout << "  LEAF took: " << timerLF.RealTime() << " sec for " << HKManager::GetME()->GetHitCollection_ID()->size() +  HKManager::GetME()->GetHitCollection_mPMT()->size() << " Hits";
			std::cout << "  (Vtx Search : "   << fLeafPerformances.vtx_search_ct 
			          << " , Vtx Minimize : " << fLeafPerformances.vtx_minimize_ct 
					  << ", Dir Search : "    << fLeafPerformances.dir_search_ct 
					  << " , Dir Minimize : " << fLeafPerformances.dir_minimize_ct 
					  << " , Energy Fit : "   << fLeafPerformances.energy_fit_ct << ")" << std::endl;
			std::cout << "  True vertex: (" << true_vertex.X()
			                        << ", " << true_vertex.Y()
									<< ", " << true_vertex.Z() 
									<< ", " << true_vertex.T()<< ")  [cm/cm/cm/ns] "
    			    << "  //  direction: (" << true_dir.X()
						            << ", " << true_dir.Y()
						            << ", " << true_dir.Z() 
					    << "  //  energy: " << true_energy << std::endl;
			std::cout << "  LEAF vertex: (" << fLeafOutput.vtx.X() 
			                        << ", " << fLeafOutput.vtx.Y()
									<< ", " << fLeafOutput.vtx.Z()
									<< ", " << fLeafOutput.vtx.T()<< ")  [cm/cm/cm/ns]" 
    			    << "  //  direction: (" << fLeafOutput.dir.X() 
						            << ", " << fLeafOutput.dir.Y()
						            << ", " << fLeafOutput.dir.Z()
					    << "  //  energy: " << fLeafOutput.energy << std::endl;
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
	fPrimaryTree->Branch("true_vertexes", &true_vertexes);
	fPrimaryTree->Branch("true_energies", &true_energies);
	fPrimaryTree->Branch("true_dir", &true_dir);

	fPrimaryTree->Branch("ID_hits", &Hit_ID, "Hit_ID/I");
	fPrimaryTree->Branch("ID_hits_50", &Hit_ID_50, "Hit_ID_50/I");
	fPrimaryTree->Branch("ID_hits_200", &Hit_ID_200, "Hit_ID_200/I");
	fPrimaryTree->Branch("ID_hits_400", &Hit_ID_400, "Hit_ID_400/I");

	fPrimaryTree->Branch("mPMT_hits", &Hit_mPMT, "Hit_mPMT/I");
	fPrimaryTree->Branch("mPMT_hits_20", &Hit_mPMT_20, "Hit_mPMT_20/I");
	fPrimaryTree->Branch("mPMT_hits_50", &Hit_mPMT_50, "Hit_mPMT_50/I");
	fPrimaryTree->Branch("mPMT_hits_200", &Hit_mPMT_200, "Hit_mPMT_200/I");
	fPrimaryTree->Branch("mPMT_hits_400", &Hit_mPMT_400, "Hit_mPMT_400/I");

	fPrimaryTree->Branch("lf_vtx", &fLeafOutput.vtx);
	fPrimaryTree->Branch("lf_vtx_nll", &fLeafOutput.vtx_nll);
	fPrimaryTree->Branch("lf_dir", &fLeafOutput.dir);
	fPrimaryTree->Branch("lf_dir_nll", &fLeafOutput.dir_nll);
	fPrimaryTree->Branch("lf_energy", &fLeafOutput.energy);
	fPrimaryTree->Branch("lf_nll", &fLeafOutput.nll_r);
	fPrimaryTree->Branch("lf_ct", &fLeafPerformances.leaf_ct);
	fPrimaryTree->Branch("lf_vtx_search_ct", &fLeafPerformances.vtx_search_ct);
	fPrimaryTree->Branch("lf_vtx_minization_ct", &fLeafPerformances.vtx_minimize_ct);
	fPrimaryTree->Branch("lf_dir_search_ct", &fLeafPerformances.dir_search_ct);
	fPrimaryTree->Branch("lf_dir_minization_ct", &fLeafPerformances.dir_minimize_ct);
	
	fPrimaryTree->Branch("rawTriggerTime", &rawTriggerTime, "rawTriggerTime/D");
	fPrimaryTree->Branch("TotalCharge", &fLeafOutput.total_charge, "TotalCharge/D");
	
	fPrimaryTree->Branch("DigiHitT", &digithit_T);
	fPrimaryTree->Branch("CorrectedDigiHitT", &correctedDigithit_T);
	fPrimaryTree->Branch("HitIsDR", &hit_is_DR);
	fPrimaryTree->Branch("DigiHitQ", &digithit_Q);
	fPrimaryTree->Branch("DigiHitResidual", &hit_residual);
	fPrimaryTree->Branch("Charge_PMT", &Charge_PMT);
}


bool AnalyseEvent(WCSimRootEvent *tEvent, PMTType pmtType) {
	const int iPMTType = static_cast<int>(pmtType);

	int nVertex = 0;				 // Number of Vertex in event
	int nTrack = 0;					 // Number of Track in event
	int nRawCherenkovHits = 0;		 // Number of Raw Cherenkov hits
	int nDigitizedCherenkovHits = 0; // Number of Digitized Cherenkov hits
	nTrigger = tEvent->GetNumberOfEvents();
	
	usedTriggerId = -1;
	rawTriggerTime = 0.;
	WCSimRootTrigger *fRootTrigger;

	//* find best trigger (best number of hits)
	int bestNumbOfHits = 0;
	bestTrigger = -1;
	for (int iTrig = 0; iTrig < nTrigger; iTrig++) {
		fRootTrigger = tEvent->GetTrigger(iTrig);
		nDigitizedCherenkovHits = fRootTrigger->GetNcherenkovdigihits();

		if(nDigitizedCherenkovHits > bestNumbOfHits && fRootTrigger->GetTriggerInfo().size() >= 3) {
			bestNumbOfHits = nDigitizedCherenkovHits;
			bestTrigger = iTrig;
		}
	}

	if(bestTrigger < 0) return false; //* Cannot process this event, none of the triggers have sufficent information

	usedTriggerId = bestTrigger;

	//* Retreive true infos (onlyh set by wcsim in the first trigger)

	for (int iTrig = 0; iTrig < tEvent->GetNumberOfEvents(); iTrig++) {
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
		if (iTrig == 0) {
			// Loop on vertex
			for (int iVertex = 0; iVertex < nVertex; iVertex++) {
#ifdef OLD_WCSIM
				true_vertex = ROOT::Math::XYZTVector(fRootTrigger->GetVtx(0),fRootTrigger->GetVtx(1),fRootTrigger->GetVtx(2),fRootTrigger->GetVtx(3));
				true_vertexes.push_back(true_vertex);
#else
				true_vertex = ROOT::Math::XYZTVector(fRootTrigger->GetVtxs(iVertex,0),fRootTrigger->GetVtxs(iVertex,1),fRootTrigger->GetVtxs(iVertex,2),fRootTrigger->GetVtxs(iVertex,3));
				true_vertexes.push_back(true_vertex);
#endif

				if (iVertex * 2 > nTrack) {
					std::cout << "ERROR: Vertex and Track number incompatible (nVertex: " << nVertex << " ; nTrack " << nTrack << ") " << std::endl;
					continue;
				}

				// Beam info are registered in the Track
				TObject *Track = (fRootTrigger->GetTracks())->At(iVertex * 2);
				WCSimRootTrack *wcTrack = dynamic_cast<WCSimRootTrack *>(Track);

				if (wcTrack->GetParentId() == 0 && wcTrack->GetIpnu()==11) {
					true_dir = ROOT::Math::XYZVector(wcTrack->GetPdir(0),wcTrack->GetPdir(1),wcTrack->GetPdir(2)).Unit();
				}
				
				true_energy = wcTrack->GetE();
				true_particleId.push_back(wcTrack->GetIpnu());
				true_energies.push_back(wcTrack->GetE());
			}

			// Send True Position to LEAF to make checks of success in the steps of leaf (should be removed for true data)

			LEAF::GetME()->SetTrueVertexInfo(true_vertexes[0]);
			LEAF::GetME()->SetTrueDirInfo(true_dir);

			rawhit_num = nRawCherenkovHits;
		}

		if(iTrig != usedTriggerId) continue; //* The trigger being processed

		auto triggerInfo = fRootTrigger->GetTriggerInfo();
		fLfTriggerTime = triggerInfo[2] - triggerInfo[1]; // time of th event - offset (common offset of all events, depends on wcsim config giles)
		rawTriggerTime = triggerInfo[2];

		startingCherenkovHitID = digithit_pmtId.size();

		// Loop on Digitized Hits
		double AllQ = 0;
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++) {
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
			digithit_Type.push_back(iPMTType);
		}


		int iIdx_BS = 0;

		nDigitizedCherenkovHits = digithit_pmtId.size();
		digithit_num = nDigitizedCherenkovHits - startingCherenkovHitID;
		std::vector<double> times;

		// Feed fitter
		if (pmtType == PMTType::kID) {
			for (int iDigitHit = startingCherenkovHitID; iDigitHit < nDigitizedCherenkovHits; iDigitHit++) {
				if (digithit_pmtId[iDigitHit] <= 0) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
				
				times.push_back(digithit_T[iDigitHit]);
				HKManager::GetME()->AddHit_ID(digithit_T[iDigitHit], digithit_Q[iDigitHit], digithit_pmtId[iDigitHit]);

				// Bonsai (do not store mPMT hits)
				if (iIdx_BS < 2000 && digithit_T[iDigitHit] < 4000. && digithit_T[iDigitHit] > -4000.) {
					// std::cout << pmtId << " " << HitT << std::endl;
					bsCAB[iIdx_BS] = digithit_pmtId[iDigitHit];
					bsT[iIdx_BS] = digithit_T[iDigitHit]; // shift BS time is needed if interaction time is < 0, this needs to be considered
					bsQ[iIdx_BS] = digithit_Q[iDigitHit];

					iIdx_BS += 1;
					bsnhit[0] = iIdx_BS;
				}
			}
		}
		else {
			for (int iDigitHit = startingCherenkovHitID; iDigitHit < nDigitizedCherenkovHits; iDigitHit++) {
				if (digithit_pmtId[iDigitHit] <= 0) std::cout << " Weird PMT ID " << digithit_pmtId[iDigitHit] << std::endl;
				times.push_back(digithit_T[iDigitHit]);
				HKManager::GetME()->AddHit_mPMT(digithit_T[iDigitHit], digithit_Q[iDigitHit], digithit_pmtId[iDigitHit]);
			}

		}

		fHit = digithit_pmtId.size();

		// Compute hits
		std::vector<double> Hit_time_20;
		std::vector<double> Hit_time_50;
		std::vector<double> Hit_time_200;
		std::vector<double> Hit_time_400;

		std::sort(times.begin(), times.end());
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++) {
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
