/*********************************************************************************/
/**	analysis.cpp								**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		**/
/**	Date: Febuary 10th 2020							**/
/**	Desc: Example application code for Benjamin's Low-E Fitter for Hyper-K	**/
/*********************************************************************************/
#include "analysis.h"


///* Some Utility Functions

double calculateToWall(double R, double h, const std::vector<double>& position, const std::vector<double>& direction) 
{
    double x = position[0], y = position[1], z = position[2];
    double dx = direction[0], dy = direction[1], dz = direction[2];

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

double calculateDWall(double R, double h, const std::vector<double>& position) 
{
    double x = position[0], y = position[1], z = position[2];
    double distanceToAxis = sqrt(x * x + y * y);
    double dSide = fabs(distanceToAxis - R);
    double dTop = fabs(z - h / 2);
    double dBottom = fabs(z + h / 2);
    return std::min({dSide, dTop, dBottom});
}

double angleBetween3DVectors(const std::vector<double>& a, const std::vector<double>& b) 
{
	if(a.size() < 3 || b.size() < 3)
	{
		throw std::invalid_argument("Vectors must be of size at least 3.");
	}
    // double dotProduct = dot(a,b);

    // Magnitudes
    double magA = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    double magB = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);

    if (magA == 0 || magB == 0) {
        throw std::invalid_argument("Zero-length vector.");
    }

    // Cosine of the angle
    double cosTheta = dot(a,b) / (magA * magB);

    // Clamp to [-1, 1] to avoid domain errors in acos
    cosTheta = std::max(-1.0, std::min(1.0, cosTheta));

    // Return angle in degrees
    return std::acos(cosTheta) * 180.0 / M_PI;
}

bool IsDarkRateHit(WCSimRootTrigger *fRootTrigger, WCSimRootCherenkovDigiHit *wcDigitHit)
{
	bool isDR = false;
	std::vector<int> rawhitphotonIDs = wcDigitHit->GetPhotonIds();
	if(rawhitphotonIDs.size() > 0)
	{
		double DR_Amount = 0;
		for(unsigned int i = 0; i < rawhitphotonIDs.size(); i++)
		{
			TObject *RawHitTimess;
			if (rawhitphotonIDs[i] >= 0 && rawhitphotonIDs[i] < (fRootTrigger->GetCherenkovHitTimes())->GetEntries()) 
			{
				RawHitTimess = (fRootTrigger->GetCherenkovHitTimes())->At(rawhitphotonIDs[i]);
			} 
			else 
			{
				std::cerr << "Error: Invalid photon ID index: " << rawhitphotonIDs[i] << std::endl;
				continue;
			}
			WCSimRootCherenkovHitTime *wcRawHitTimee = dynamic_cast<WCSimRootCherenkovHitTime *>(RawHitTimess);
			// std::cout << "Photon type : " << wcRawHitTimee->GetPhotonCreatorProcessName() << std::endl;
			if(wcRawHitTimee->GetPhotonCreatorProcessName() == "darkNoise") DR_Amount+=1.0;
		}
		DR_Amount = DR_Amount / rawhitphotonIDs.size();
		if(DR_Amount > 0.1) isDR = true; //? what threshold to consider a hit as darkRate ?
	}
	else
	{
		std::cout << "Error: No photon IDs found for this hit : " << wcDigitHit->GetTubeId() << std::endl;
	}
	return isDR;
}

double GetResidualTime(const std::vector<double>& origin, double originTime, WCSimRootCherenkovDigiHit *wcDigitHit, TimeDelta fTimeCorrection, double fLfTriggerTime)
{
	float fCVacuum = 3e8 * 1e2 / 1e9; // speed of light, in centimeter per ns.
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

double GetDistanceToNeighbors(WCSimRootTrigger *fRootTrigger, WCSimRootCherenkovDigiHit *wcDigitHit, int N_Neighbors)
{
	std::vector<double> distances;

	WCSimRootPMT pmt;
	pmt = fLeafGeometry->GetPMT(wcDigitHit->GetTubeId() - 1, false);

	std::vector<double> PMTpos(3);
	for (int j = 0; j < 3; j++) PMTpos[j] = pmt.GetPosition(j);

	//loop over each hit, check if it is not the same as the one we are looking at
	for (int i = 0; i < fRootTrigger->GetNcherenkovdigihits(); i++)
	{
		WCSimRootCherenkovDigiHit *neighborHit =  dynamic_cast<WCSimRootCherenkovDigiHit *>((fRootTrigger->GetCherenkovDigiHits())->At(i));
			
		if(wcDigitHit->GetTubeId() == neighborHit->GetTubeId())
			continue;

		WCSimRootPMT neighborPMT;
		neighborPMT = fLeafGeometry->GetPMT(neighborHit->GetTubeId() - 1, false);

		std::vector<double> neighborPMTpos(3);
		for (int j = 0; j < 3; j++) neighborPMTpos[j] = neighborPMT.GetPosition(j);

		double distance = calculateDistance(PMTpos, neighborPMTpos);
		distances.push_back(distance);
	}

	//* Sort the distances, keep the N_Neighbors closest ones
	std::sort(distances.begin(), distances.end());
	//  size = N_Neighbors;
	if (distances.size() < (long unsigned int)N_Neighbors )
	{
		std::cerr << "Not enough neighbors found!" << std::endl;
		return -1;
	}
	double sum = 0;
	for (int i = 0; i < N_Neighbors; i++) sum += distances[i];
	return  sum / N_Neighbors;
}

///* MAIN EXECUTION

int main(int argc, char **argv)
{
	std::string sInputFile = "";
	std::string sOutputFile = "";

	int iNeededArgc = 3;
	double dDarkNoise = 0.;		  // kHz
	double dDarkNoiseHybrid = 0.; // kHz

	iNeededArgc += 1;
#ifdef mPMT
	iNeededArgc += 1;
#endif

	if (argc == iNeededArgc)
	{
		int iArg = 1;
		sInputFile = argv[iArg];
		iArg += 1;
		sOutputFile = argv[iArg];
		iArg += 1;
		dDarkNoise = atof(argv[iArg]);
		iArg += 1;
#ifdef mPMT
		dDarkNoiseHybrid = atof(argv[iArg]);
		iArg += 1;
#endif
	}
	else
	{

		std::cout << "Synthax: " << argv[0] << " input output";
		std::cout << " DN_in_kHz_B&L";

#ifdef mPMT
		std::cout << " DN_in_kHz_mPMT";
#endif
		std::cout << std::endl;
		return 0;
	}

	// Read WCSim output
	TFile *fInputFile = new TFile(sInputFile.c_str(), "READ");

	// Get TTrees
	TTree *fInputTree = (TTree *)fInputFile->Get("wcsimT");
	TTree *fInputGeoTree = (TTree *)fInputFile->Get("wcsimGeoT");

	fLeafGeometry = 0;
	fInputGeoTree->SetBranchAddress("wcsimrootgeom", &fLeafGeometry);

	WCSimRootEvent *fIDevent = new WCSimRootEvent();
	// Set Branch
	fInputTree->SetBranchAddress("wcsimrootevent", &fIDevent);
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	// WCSimRootEvent *fHybridevent = new WCSimRootEvent();
#ifdef mPMT
	// Set Branche
	// fInputTree->SetBranchAddress("wcsimrootevent2", &fHybridevent);
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
#endif

#ifdef OD_ON
	WCSimRootEvent *fODevent = new WCSimRootEvent();
	// Set Branche
	fInputTree->SetBranchAddress("wcsimrootevent_OD", &fODevent);
	// Set autodelete to avoid memory leak
	fInputTree->GetBranch("wcsimrootevent_OD")->SetAutoDelete(kTRUE);
#endif

	// Read Geo
	fInputGeoTree->GetEntry(0);

	// Create Output TTree
	TFile *fOutputFile = new TFile(sOutputFile.c_str(), "RECREATE");
	fOutputFile->SetCompressionLevel(2);
	TTree *fGeoTree = new TTree("wcsimGeoT", "Geometry TTree");
	TTree *fPrimaryTree = new TTree("Reduced", "Reduced TTree");

	// Set Branches
	SetGeoBranch(fGeoTree);

	FitterOutput leaf_output;
	SetCustomBranch(fPrimaryTree, leaf_output);

	fGeoTree->Fill();

	// Get PMT Number:
	int nPMT_ID = fLeafGeometry->GetWCNumPMT();
#ifdef OD_ON
	int nPMT_OD = fLeafGeometry->GetODWCNumPMT();
#endif
	int nMultPMT = fLeafGeometry->GetWCNumPMT(true);

	std::cout << " ID " << nPMT_ID << std::endl;
	std::cout << " mPMT " << nMultPMT << std::endl;

	HKManager::GetME()->SetGeometry(fLeafGeometry, dDarkNoise * 1e3, dDarkNoiseHybrid * 1e3);

	// Initialize LEAF
	LEAF::GetME()->Initialize(HKManager::GetME()->GetGeometry()); // Loads the geometry in leaf and loads the pdfs

	//* for each config variable, check if it has been set, take the default variable if not
	int maxEvents = -1;		 // in cm, the step size for coarse grid search
	const char* envValue = std::getenv("nbOfEvents");
	if (envValue != nullptr) 
	{
		maxEvents = std::stoi(envValue);
		std::cout << "max events is set to " << maxEvents << std::endl;
	} else std::cout << "Environment variable max events is not set. Using all events" << std::endl;

	envValue = std::getenv("maxHitsAngle");
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
	int nPrimaryEvents = fInputTree->GetEntries();

	if(maxEvents > 0) nPrimaryEvents = std::min(nPrimaryEvents, maxEvents);

	int iWrite = 0;
	int failAmount = 0;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	// Loop on Primary events
	for (int i = 0; i < nPrimaryEvents; i++)
	{
		// Reset Hit vector
		HKManager::GetME()->ResetHitInfo();
		// HKManager::GetME()->ResetSecondaryHitInfo();

		// if ( i%1000==0 ) {
		timer.Stop();
		std::cout << "Event # = " << i << " / " << nPrimaryEvents << " ( " << timer.RealTime() << " )\n";
		timer.Reset();
		timer.Start();
		//}

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
		digithit_Angle.clear();
		digithit_NormAngle.clear();
		digithit_NeighborsDist.clear();
		hit_is_DR.clear();
		Charge_PMT.clear();
		relativeAngle.clear();
		lf_relativeAngle.clear();
		hit_residual.clear();

		rawhit_num = 0;
		digithit_num = 0;

		fLFTime = -9999;

		hit_5_15ns = 0;
		hit_50ns = 0;

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

		lf_spatial_res = 0.;
		lf_Dir_res = 0.;
		lf_time_res = 0.;
		lf_energy_res = 0.;

		bestTrigger = 0;
		fLfTriggerTime = 0.;

		/****************************************************************************************/
		/* ID events										*/
		/****************************************************************************************/

		fHit = 0;
		fHit_20 = 0;
		fHit_50 = 0;
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

		fHit = 0;
		fHit_20 = 0;
		fHit_50 = 0;
		fHit_200 = 0;
		fHit_400 = 0;

#ifdef mPMT
		// /*bool bmPMT =*/AnalyseEvent(fHybridevent, mPMT_EVENT);
#endif
		Hit_mPMT = fHit;
		Hit_mPMT_20 = fHit_20;
		Hit_mPMT_50 = fHit_50;
		Hit_mPMT_200 = fHit_200;
		Hit_mPMT_400 = fHit_400;

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
		/*bool bOD =*/AnalyseODEvent(fODevent, OD_EVENT);

		Hit_OD = fHit;
		Hit_OD_50 = fHit_50;
		Hit_OD_200 = fHit_200;
		Hit_OD_400 = fHit_400;
#endif
		/****************************************************************************************/
		/* Benjamin Fitter									*/
		/****************************************************************************************/

		TStopwatch timerLF;
		timerLF.Reset();
		timerLF.Start();

		// To be replaced by meaningful TriggerTime
		TimeDelta fDummyTrigger(0.);

		leaf_output = LEAF::GetME()->MakeSequentialFit(HKManager::GetME()->GetHitCollection(), fDummyTrigger);
		// the fit method can be replaced by a joint fit (doesn't work better for now)
		// leaf_output = LEAF::GetME()->MakeJointFit(HKManager::GetME()->GetHitCollection(), fDummyTrigger);
		fOutputProps = LEAF::fOutputProps;
		
		timerLF.Stop();

		fLFTime = timerLF.RealTime();

		std::cout << " LEAF took: " << timerLF.RealTime() << " for " << HKManager::GetME()->GetHitCollection()->Size() << " Hits";
		std::cout << " (Vtx Search : " << fOutputProps.Vtx_Search_ComputeTime << " , Vtx Minimize : " << fOutputProps.Vtx_Minimize_ComputeTime << ", Dir Search : " << fOutputProps.Dir_Search_ComputeTime << " , Dir Minimize : " << fOutputProps.Dir_Minimize_ComputeTime << " , Energy Fit : " << fOutputProps.Energy_Fit_ComputeTime << ")" << std::endl;
		bool validEvent = PostLeafAnalysis(fIDevent,ID_EVENT, leaf_output);

		/****************************************************************************************/
		/* Fill output tree									*/
		/****************************************************************************************/
		if(validEvent)
		{
			fPrimaryTree->Fill();
			iWrite += 1;
		}
	}

	std::cout << "Event # = " << nPrimaryEvents << " / " << nPrimaryEvents << std::endl;

	std::cout << ((nPrimaryEvents - failAmount) / nPrimaryEvents) * 100 << "% of the events were processed" << std::endl;

	fOutputFile->cd(); 

	fOutputFile->Write("", TObject::kOverwrite);

	delete fPrimaryTree;
	delete fOutputFile;
	return 1;
}

//* All the variables that will be kept in the output Tree
void SetCustomBranch(TTree *fPrimaryTree, FitterOutput leaf_output)
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
	fPrimaryTree->Branch("DWall", &dWall, "DWall/D");
	fPrimaryTree->Branch("ToWall", &toWall, "lf_wall/D");

	fPrimaryTree->Branch("ID_hits", &Hit_ID, "Hit_ID/I");
	fPrimaryTree->Branch("ID_hits_50", &Hit_ID_50, "Hit_ID_50/I");
	fPrimaryTree->Branch("ID_hits_200", &Hit_ID_200, "Hit_ID_200/I");
	fPrimaryTree->Branch("ID_hits_400", &Hit_ID_400, "Hit_ID_400/I");

	fPrimaryTree->Branch("mPMT_hits", &Hit_mPMT, "Hit_mPMT/I");
	fPrimaryTree->Branch("mPMT_hits_20", &Hit_mPMT_20, "Hit_mPMT_20/I");
	fPrimaryTree->Branch("mPMT_hits_50", &Hit_mPMT_50, "Hit_mPMT_50/I");
	fPrimaryTree->Branch("mPMT_hits_200", &Hit_mPMT_200, "Hit_mPMT_200/I");
	fPrimaryTree->Branch("mPMT_hits_400", &Hit_mPMT_400, "Hit_mPMT_400/I");

	fPrimaryTree->Branch("hit_5_15ns", &hit_5_15ns, "hit_5_15ns/I");
	fPrimaryTree->Branch("hit_50ns", &hit_50ns, "hit_50ns/I");

	fPrimaryTree->Branch("lf_vertex", &leaf_Vertex);
	fPrimaryTree->Branch("lf_NLL", &leaf_output.NLL, "lf_NLL/D");
	fPrimaryTree->Branch("lf_good", &leaf_output.NLLR, "lf_good/D");
	fPrimaryTree->Branch("lf_ctime", &fLFTime, "lf_ctime/D"); // Computation time
	fPrimaryTree->Branch("lf_energy", &leaf_output.Energy, "lf_energy/D");
	fPrimaryTree->Branch("lf_Dir", &leaf_Dir);
	fPrimaryTree->Branch("lf_MyDir", &leaf_MyDir);
	fPrimaryTree->Branch("lf_Quick_Dir", &leaf_QuickDir);
	fPrimaryTree->Branch("lf_Dir_NLL", &leaf_output.DNLL, "lf_Dir_NLL/D");
	fPrimaryTree->Branch("lf_DWall", &lf_dWall, "DWall/D");
	fPrimaryTree->Branch("lf_ToWall", &lf_ToWall, "lf_ToWall/D");

	fPrimaryTree->Branch("lf_spatial_res", &lf_spatial_res, "lf_spatial_res/D");
	fPrimaryTree->Branch("lf_time_res", &lf_time_res, "lf_time_res/D");
	fPrimaryTree->Branch("lf_Dir_res", &lf_Dir_res,"lf_Dir_res/D");
	fPrimaryTree->Branch("lf_MyDir_res", &lf_MyDir_res,"lf_MyDir_res/D");
	fPrimaryTree->Branch("lf_Quick_Dir_res", &lf_Quick_Dir_res,"lf_Quick_Dir_res/D");
	fPrimaryTree->Branch("lf_energy_res", &lf_energy_res, "lf_energy_res/D");

	fPrimaryTree->Branch("Leaf_ComputeTime", &fOutputProps.Leaf_ComputeTime, "Leaf_ComputeTime/D");
	fPrimaryTree->Branch("Vtx_Search_ComputeTime", &fOutputProps.Vtx_Search_ComputeTime, "Vtx_Search_ComputeTime/D");
	fPrimaryTree->Branch("Vtx_Minimize_ComputeTime", &fOutputProps.Vtx_Minimize_ComputeTime, "Vtx_Minimize_ComputeTime/D");
	fPrimaryTree->Branch("Dir_Quick_Search_ComputeTime", &fOutputProps.Dir_Quick_Search_ComputeTime, "Dir_Quick_Search_ComputeTime/D");
	fPrimaryTree->Branch("Dir_Search_ComputeTime", &fOutputProps.Dir_Search_ComputeTime, "Dir_Search_ComputeTime/D");
	fPrimaryTree->Branch("Dir_Minimize_ComputeTime", &fOutputProps.Dir_Minimize_ComputeTime, "Dir_Minimize_ComputeTime/D");
	fPrimaryTree->Branch("Energy_Fit_ComputeTime", &fOutputProps.Energy_Fit_ComputeTime, "Energy_Fit_ComputeTime/D");

	fPrimaryTree->Branch("rawTriggerTime", &rawTriggerTime, "rawTriggerTime/D");
	fPrimaryTree->Branch("TotalCharge", &leaf_output.TotalCharge, "TotalCharge/D");
	fPrimaryTree->Branch("RelativeAngle", &relativeAngle);
	fPrimaryTree->Branch("lf_RelativeAngle", &lf_relativeAngle);
	
	fPrimaryTree->Branch("DigiHitT", &digithit_T);
	fPrimaryTree->Branch("CorrectedDigiHitT", &correctedDigithit_T);
	fPrimaryTree->Branch("HitIsDR", &hit_is_DR);
	fPrimaryTree->Branch("DigiHitQ", &digithit_Q);
	fPrimaryTree->Branch("DigiHitAngle", &digithit_Angle);
	fPrimaryTree->Branch("DigiHitNormAngle", &digithit_NormAngle);
	fPrimaryTree->Branch("DigiHitNeighborsDist", &digithit_NeighborsDist);
	fPrimaryTree->Branch("DigiHitResidual", &hit_residual);
	fPrimaryTree->Branch("Charge_PMT", &Charge_PMT);
}

void SetGeoBranch(TTree *fGeoTree)
{
	fGeoTree->Branch("wcsimrootgeom", fLeafGeometry);
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

				toWall = calculateToWall(R, h, trueVertex, trueDir);
				dWall = calculateDWall(R, h, trueVertex);

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

		//* check if trigger is too soon or too late compared to the event (trigger happened because of dark rate)
		double hitTimesStack = 0;
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++)
		{
			TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
			WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit *>(Hit);
			double HitT = wcDigitHit->GetT() + fLfTriggerTime;
			hitTimesStack+=HitT;
		}
		double hitTmean = hitTimesStack/nDigitizedCherenkovHits;
		if(hitTmean > 2000.0 || hitTmean < -2000.0) return false; //* set these values to define the filter (currently filtering nothing)

		double AllQ = 0;

		// Loop on Digitized Hits
		for (int iDigitHit = 0; iDigitHit < nDigitizedCherenkovHits; iDigitHit++)
		{
			//* Get and Compute Hit Informations

			TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
			WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit *>(Hit);
			
			int pmtId = wcDigitHit->GetTubeId();
			WCSimRootPMT pmt = fLeafGeometry->GetPMT(pmtId - 1, false);
			double HitT = wcDigitHit->GetT() + fLfTriggerTime;
			double HitQ = wcDigitHit->GetQ();
			int peForTube = wcDigitHit->GetQ();
			std::vector<double>  lfHitDir = std::vector<double>(3);
			std::vector<double> PMTOrientation = std::vector<double>(3);
			std::vector<double> vDir = std::vector<double>(3);
			bool isDR = IsDarkRateHit(fRootTrigger, wcDigitHit);
			double hitResidual = GetResidualTime(trueVertex, true_origin_T[0], wcDigitHit, fTimeCorrection, fLfTriggerTime);
			for (int j = 0; j < 3; j++) 
			{
				PMTOrientation[j] = pmt.GetOrientation(j);
				vDir[j] = pmt.GetPosition(j) - trueVertex[j];
				lfHitDir[j] = pmt.GetPosition(j);
			}
			
			Normalize(vDir);
			Normalize(lfHitDir);
			Normalize(PMTOrientation);

			std::vector<double> neg_vDir = { -vDir[0], -vDir[1], -vDir[2] };
			double NormAngle = dot(neg_vDir, PMTOrientation);

			double hitRelativeAngle = dot(vDir, trueDir);

			//* Use that inforation to filter some hits
			//* we aim at killing dark rate hits
			//* ideally without using any information form the true event
			
			//* Direction Filter
			double angle = angleBetween3DVectors(vDir, trueDir);
			// if(angle > maxHitAngle) continue; // works best at 90°

			//* Filter base on distance to hit neighbors
			double distanceToNeighbor = GetDistanceToNeighbors(fRootTrigger, wcDigitHit, N_Neighbors);
			// if(distanceToNeighbor > maxDistanceToNeighbors) continue;
			
			//*DR FILTER
			if(isDR) continue; // ideal filter to check reconstruction without dark rate
			
			///* fill infos to output root file
			AllQ += HitQ;
			digithit_pmtId.push_back(pmtId);
			digithit_T.push_back(HitT);
			correctedDigithit_T.push_back(HitT - fLfTriggerTime - triggerInfo[1]);
			digithit_Q.push_back(HitQ);
			hit_is_DR.push_back(isDR);
			digithit_Angle.push_back(angle);
			digithit_NormAngle.push_back((TMath::ACos(NormAngle))*180./TMath::Pi());
			digithit_NeighborsDist.push_back(distanceToNeighbor);
			Charge_PMT.push_back(peForTube);
			relativeAngle.push_back((TMath::ACos(hitRelativeAngle))*180./TMath::Pi());
			digithit_Type.push_back(iEventType);
			hit_residual.push_back(hitResidual);
		}

		for(long unsigned int j=0;j<hit_residual.size();j++)
		{
 			if(hit_residual[j] < 15 && hit_residual[j] > -5) hit_5_15ns++;
			if(hit_residual[j] < 50 && hit_residual[j] > -50) hit_50ns++;
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

			// Count hit in 200 ns windowAnalyseEvent
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

bool PostLeafAnalysis(WCSimRootEvent * tEvent, int iEventType, FitterOutput leaf_output)
{
	if (!(true_origin_X.size() > 0 && true_origin_Y.size() > 0 && true_origin_Z.size() > 0)) return false;

	for(int j=0;j<4;j++) leaf_Vertex.push_back(leaf_output.Vtx[j]);
	for(int j=0;j<3;j++) leaf_Dir.push_back(leaf_output.Dir[j]);
	for(int j=0;j<3;j++) leaf_MyDir.push_back(leaf_output.MyDir[j]);
	for(int j=0;j<3;j++) leaf_QuickDir.push_back(leaf_output.Quick_Dir[j]);

	//* RESOLUTIONS COMPUTATION
	std::vector<double> true_vertex = {true_origin_X[0], true_origin_Y[0], true_origin_Z[0]};
	lf_dWall = calculateDWall(R, h, leaf_output.Vtx);
	lf_ToWall = calculateToWall(R, h,leaf_output.Vtx, leaf_output.Dir);
	lf_spatial_res = calculateDistance(true_vertex, leaf_output.Vtx);
	lf_time_res = abs(leaf_output.Vtx[3]-true_origin_T[0]);

	if(leaf_output.Dir.size() == 3) lf_Dir_res = TMath::ACos(dot(leaf_output.Dir, trueDir))*180./TMath::Pi();
	if(leaf_output.MyDir.size() == 3) lf_MyDir_res = TMath::ACos(dot(leaf_output.MyDir, trueDir))*180./TMath::Pi();
	if(leaf_output.Quick_Dir.size() == 3) lf_Quick_Dir_res = TMath::ACos(dot(leaf_output.Quick_Dir, trueDir))*180./TMath::Pi();

	lf_energy_res = abs(leaf_output.Energy - true_energy[0]);

	//* ANGLE ANALYSIS
	WCSimRootTrigger * fRootTrigger = tEvent->GetTrigger(bestTrigger);

	for(int iDigitHit = startingCherenkovHitID; (long unsigned int)iDigitHit < digithit_pmtId.size(); iDigitHit++)
	{
		TObject *Hit = (fRootTrigger->GetCherenkovDigiHits())->At(iDigitHit);
		WCSimRootCherenkovDigiHit *wcDigitHit = dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
		WCSimRootPMT pmt = fLeafGeometry->GetPMT(wcDigitHit->GetTubeId() - 1, false);

		std::vector<double> lfHitDir = std::vector<double>(3);
		for(int j=0;j<3;j++) lfHitDir[j] = pmt.GetPosition(j) - leaf_output.Vtx[j];
		Normalize(lfHitDir);

		std::vector<double> trueHitDir = std::vector<double>(3);
		for(int j=0;j<3;j++) trueHitDir[j] = pmt.GetPosition(j) - trueVertex[j];
		Normalize(trueHitDir);
		double hitlfRelativeAngle = dot(trueHitDir, leaf_output.Dir);
		lf_relativeAngle.push_back((TMath::ACos(hitlfRelativeAngle))*180./TMath::Pi());
	}
	return true;
}