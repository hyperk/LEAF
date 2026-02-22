/*****************************************************************************************************/
/**	LEAF.cc											**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)					**/
/**	Original author: Benjamin Quilain								**/
/**	Date: December 18th 2019									**/
/**	Desc: Low-E Fitter for Hyper-K								**/
/*****************************************************************************************************/

#include "LEAF.hh"
#include "TStopwatch.h"

LEAF *LEAF::myFitter = NULL;
std::mutex mtx;
FitterOutput LEAF::fOutput;
FitterOutputProps LEAF::fOutputProps;


/*****************************************************************************************************/
/* GENERAL LEAF METHODS AND STRUCTURES FOR FITS */
/*****************************************************************************************************/

LEAF::LEAF()
{
	myFitter = this;
}

LEAF::~LEAF()
{
	delete fRand;
}

LEAF *LEAF::GetME()
{
	if (myFitter)
		return myFitter;

	myFitter = new LEAF();
	return myFitter;
}

void LEAF::DeleteME()
{
	if (myFitter) delete myFitter;
}

void LEAF::Initialize(const Geometry *lGeometry)
{
	InitConfig(lGeometry);
	InitInputs(lGeometry);
	LoadSplines();
	fOutput = NewOutuput();
}

void LEAF::LoadHitCollection(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT)
{
	fHitCollection = lHitCol;
	fTriggerTime = lTriggerTime;
	fTimeCorrection = fHitCollection->timestamp - fTriggerTime;
	if (fPositionList.size() == 0) MakePositionList();
	fUseDirectionality = bMultiPMT ? fUseDirectionality : false;
}

struct FitterOutput LEAF::NewOutuput()
{
	struct FitterOutput fOutput;
	fOutput.Vtx = std::vector<double>(4);
	fOutput.Dir = std::vector<double>(3);
	fOutput.NLL		= 0.;
	fOutput.DNLL	= 0.;
	fOutput.Energy = 0.;
	fOutput.TotalCharge = 0.;

	return fOutput;
}

struct FitterOutput LEAF::MakeSequentialFit(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT)
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	LoadHitCollection(lHitCol, lTriggerTime, bMultiPMT);
	
	if(fHitCollection->Size() <= 0) return NewOutuput();
	
	//* Vertex
	FitVertex();

	//* Goodness of fit
	fOutput.NLLR = Likelihoods::GoodnessOfFit(fHitCollection, fOutput.Vtx, fHitCollection->Size(), fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality); //* Goodness of Fit
	
	//* Direction

	//? three direction methods are stored, one can use only the Fit Direction method
	// FitDirectionQuick(fOutput.Vtx); // prefit

	// // the prefit is stored in fOutput.Dir
	// FitDirection(fOutput.Vtx, fHitCollection->Size(), true); // optimization with one candidate that is the prefit
	// fOutput.MyDir = fOutput.Dir;
	// FitDirection(fOutput.Vtx, fHitCollection->Size(), false); // optimization without prefit, 30 candidates used

	FitDirectionQuick(fTrueVtxPos); // prefit

	// the prefit is stored in fOutput.Dir
	FitDirection(fTrueVtxPos, fHitCollection->Size(), true); // optimization with one candidate that is the prefit
	fOutput.MyDir = fOutput.Dir;
	FitDirection(fTrueVtxPos, fHitCollection->Size(), false); // optimization without prefit, 30 candidates used
	

	//* Energy
	FitEnergy();

	timer.Stop();
	fOutputProps.Leaf_ComputeTime = timer.RealTime();

	return fOutput;
}

struct FitterOutput LEAF::MakeJointFit(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT)
{
	LoadHitCollection(lHitCol, lTriggerTime, bMultiPMT);

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	if(fHitCollection->Size() <= 0) return NewOutuput();

	std::vector<JointFitCandidate> Candidates = SearchVertexAndDir();
	JointFitCandidate finalCandidate = MinimizeVertexAndDir(Candidates);

	fOutput.Vtx  = {finalCandidate.VtxPart.X, finalCandidate.VtxPart.Y, finalCandidate.VtxPart.Z, finalCandidate.VtxPart.T};
	fOutput.Dir  = PolarToCartesianNorm({finalCandidate.DirPart.theta, finalCandidate.DirPart.phi});
	fOutput.NLL  = finalCandidate.NLL;
	fOutput.DNLL = Likelihoods::Dir_NLL(fHitCollection, fOutput.Vtx, finalCandidate.DirPart.theta, finalCandidate.DirPart.phi , fHitCollection->Size());
	fOutput.NLLR = Likelihoods::GoodnessOfFit(fHitCollection, fOutput.Vtx, fHitCollection->Size(), fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality); //* Goodness of Fit
	
	//* Energy, alone
	FitEnergy();

	timer.Stop();
	fOutputProps.Leaf_ComputeTime = timer.RealTime();

	return fOutput;
}


/*****************************************************************************************************/
/* DIRECTION FIT */
/*****************************************************************************************************/

//Creates the candidates list after Coarse grid search, refines the fit, then outputs the best candidate
std::vector<double> LEAF::FitDirection(const std::vector<double>& fixedVertexPosition, int nhits, bool searchPrior) 
{
	//* Init
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	double *tLimits = new double[4];
	int verbose = VERBOSE;
	for (int i = 0; i < 4; i++) tLimits[i] = 2 * fSearchVtxStep;

	//* Fit first step : find candidates for the direction
	std::vector<DirectionCandidate> candidates;
	if(searchPrior)
	{
		//* Only one candidate, chosen wisely by computing a mean hit direction for hits within a specified residual time window
		std::vector<double> dirPriorCart = FitDirectionQuick(fixedVertexPosition);
		std::vector<double> dirPrior = CartesianToPolarNorm(dirPriorCart);
		DirectionCandidate priorCandidate;
		priorCandidate.theta = dirPrior[0];
		priorCandidate.phi = dirPrior[1];
		candidates.push_back(priorCandidate);
	}
	else candidates = SearchDirection(fixedVertexPosition, DirTolerance); //* multiple candidates, spread uniformly

	// std::cout << "nb direction candidates : " << candidates.size() << std::endl;

	fOutputProps.Dir_Search_ComputeTime = timer.RealTime();
	TStopwatch timer2;
	timer2.Reset();
	timer2.Start();

	//? To be removed for real data where true dir is not known, no impact on the fit
	// fOutput.stepOneContainsTrueDir = (int)ContainsTrueDir(&candidates, fTrueDir); //* a quick check to see if the step one works (one of the candidate is near the true direction)

    if (verbose >= 2) 
	{
        for (unsigned int j = 0; j < candidates.size(); j++) std::cout << "Candidate theta = " << candidates[j].theta * TMath::RadToDeg() << " deg, phi = " << candidates[j].phi * TMath::RadToDeg() << " deg, NLL = " << candidates[j].DNLL << std::endl;
    }

	DirectionCandidate finalCandidate = MinimizeDirection(fixedVertexPosition, candidates, nhits, tLimits, verbose);

	fOutput.Dir = PolarToCartesianNorm({finalCandidate.theta, finalCandidate.phi});
	fOutput.DNLL = finalCandidate.DNLL;

	fOutputProps.Dir_Minimize_ComputeTime = timer2.RealTime();

    return fOutput.Dir;
}

std::vector<double> LEAF::FitDirectionQuick(const std::vector<double>& fixedVertexPosition)
{
	// Quick fit implementation

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	const std::vector<double>& UsedVertex = fixedVertexPosition;

	std::vector<std::vector<double>> toHitVectors;
	std::vector<double> residuals;

	for(unsigned int i = 0; i < fHitCollection->Size(); i++)
	{
		std::vector<double> pmtPosition;
		std::vector<double> toHit;
		
		Hit lHit = fHitCollection->At(i);
		
		int iPMT = lHit.PMT;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];

		for(int j = 0; j < 3 ; j++) pmtPosition.push_back(lPMTInfo.Position[j]);
		for(int j = 0; j < 3 ; j++) toHit.push_back(pmtPosition[j] - UsedVertex[j]);
		Normalize(toHit);
		for(int j = 0; j < 3 ; j++) toHit[j] *= lHit.Q; // Scale by charge
		toHitVectors.push_back(toHit);

		residuals.push_back(ComputeResidualTime(UsedVertex, UsedVertex[3], lHit));
	}
	std::vector<double> lf_Dir(3, 0.0);
	for(long unsigned int i = 0; i < toHitVectors.size(); i++)
	{
		if(residuals[i] < fDirectionPDF_maxResidual && residuals[i] > fDirectionPDF_minResidual) for(int j = 0; j < 3; j++) lf_Dir[j] += toHitVectors[i][j];
	}

	fOutput.Dir = Normalize(lf_Dir);
	fOutput.Quick_Dir = Normalize(lf_Dir);
	fOutputProps.Dir_Quick_Search_ComputeTime = timer.RealTime();
	return fOutput.Dir;
}

std::vector<DirectionCandidate> LEAF::SearchDirection(const std::vector<double>& fixedVertexPosition, int tolerance)
{
	std::vector<DirectionCandidate> candidates;

    // Grid search over theta and phi
    for (double theta = 0.0; theta <= TMath::Pi(); theta += theta_step) {
        for (double phi = -TMath::Pi(); phi < TMath::Pi(); phi += phi_step) {
            // Compute NLL for this direction
            double DNLL = Likelihoods::Dir_NLL(fHitCollection ,fixedVertexPosition, theta, phi, fHitCollection->Size());

            // Store the candidate
            DirectionCandidate candidate;
            candidate.theta = theta;
            candidate.phi = phi;
            candidate.DNLL = DNLL;
            candidates.push_back(candidate);
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const DirectionCandidate& a, const DirectionCandidate& b) 
	{
        return a.DNLL < b.DNLL;
    });

	// Select the best candidates based on tolerance
	while((int)candidates.size() > tolerance) candidates.pop_back();

	return candidates;
}

//MIGRAD Optimization
DirectionCandidate LEAF::MinimizeDirection(const std::vector<double>& fixedVertexPosition, std::vector<DirectionCandidate>& candidates, int nhits,  double *limits, int verbose) 
{
    // Create a minimizer
    TFitter* minimizer = new TFitter(7); // 2 parameters: theta, phi
    TMinuit* minuit = minimizer->GetMinuit();

    double arglist[10];
    int err = 0;
    double p1 = verbose - 1;
    minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);
    minuit->SetErrorDef(1);

    if (verbose < 2) minuit->mnexcm("SET NOWarnings", 0, 0, err);

    arglist[0] = 2;
    minuit->mnexcm("SET STR", arglist, 1, err);
    minimizer->SetFCN(MinuitDirNLL);

	for (int icand = 0; icand < (int)candidates.size(); icand++)
	{
		// Set initial parameters and bounds
		minimizer->SetParameter(0, "theta", candidates[icand].theta, 2*M_PI/180, 0.0, TMath::Pi());
		minimizer->SetParameter(1, "phi", candidates[icand].phi, 2*M_PI/180, -TMath::Pi(), TMath::Pi());
		minimizer->SetParameter(2, "nhits", nhits, nhits, nhits - 1, nhits + 1);
		minimizer->SetParameter(3, "vertex0", fixedVertexPosition[0], 5e-1, -limits[0], limits[0]);
		minimizer->SetParameter(4, "vertex1", fixedVertexPosition[1], 5e-1, -limits[1], limits[1]);
		minimizer->SetParameter(5, "vertex2", fixedVertexPosition[2], 5e-1, -limits[2], limits[2]);
		minimizer->SetParameter(6, "vertex3", fixedVertexPosition[3], 5e-1 / fLightSpeed, -limits[3] / fLightSpeed, limits[3] / fLightSpeed);
		
		minimizer->FixParameter(2);
		minimizer->FixParameter(3);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);

		// Minimize
		arglist[0] = 1e9; // Maximum number of function evaluations
		arglist[1] = 1e-1; // Tolerance
		minimizer->ExecuteCommand("MIGRAD", arglist, 2 );

		// Retrieve results
		std::vector<double> optimizedParameters(2);
		candidates[icand].theta = minimizer->GetParameter(0); // theta
		candidates[icand].phi = minimizer->GetParameter(1); // phi

		minuit = minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);

		candidates[icand].DNLL = fmin;
	}

	std::sort(candidates.begin(), candidates.end(), [](const DirectionCandidate& a, const DirectionCandidate& b) 
	{
        return a.DNLL < b.DNLL;
    });
    // Optional: Print the refined parameters
    // if (verbose >= 2) std::cout << "After refinement: theta = " << optimizedParameters[0] * TMath::RadToDeg() << " deg, phi = " << optimizedParameters[1] * TMath::RadToDeg() << " deg" << std::endl;

    delete minimizer;
    return candidates[0]; // Returns [theta, phi]
}


/*****************************************************************************************************/
/* VERTEX FIT */
/*****************************************************************************************************/


void LEAF::FitVertex(std::vector<double>* fDirection_Filter, float FilterThreshold)
{
	double *tLimits = new double[4];
	int iHitsTotal = fHitCollection->Size();
	std::vector< std::vector<double> > fRecoVtxPosFinal;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<VtxCandidate> tRecoVtxPosCand = this->SearchVertex_Main(iHitsTotal,fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
	std::vector<std::vector<double>> tRecoVtxPos;

	// Extract only the vertex positions from the VtxCandidate structure
	for (const auto& cand : tRecoVtxPosCand) 
	{
		tRecoVtxPos.push_back({cand.X, cand.Y, cand.Z, cand.T, cand.NLL});
	}

	if (VERBOSE >= 2)
	{
		for (unsigned int a = 0; a < tRecoVtxPos.size(); a++)
			std::cout << "Candidate vertex [ " << a << " ] = " << tRecoVtxPos[a][0] << " , " << tRecoVtxPos[a][1] << " , " << tRecoVtxPos[a][2] << " , " << tRecoVtxPos[a][3] << ", NLL = " << tRecoVtxPos[a][4] << std::endl;
	}

	timer.Stop();
	fOutputProps.Vtx_Search_ComputeTime = timer.RealTime();
	// std::cout << "SearchVertex took: " << timer.RealTime() << std::endl;

	double dStepSizeFinal = 10;
	int iToleranceFinal = 1;

	timer.Reset();
	timer.Start();
	
	for (int i = 0; i < 4; i++) tLimits[i] = 2 * fSearchVtxStep;
	fRecoVtxPosFinal = this->MinimizeVertex_Main(tRecoVtxPos, tLimits, dStepSizeFinal, iHitsTotal, fSearchVtxTolerance, iToleranceFinal, VERBOSE, true, false, fMinimizeLimitsNegative, fMinimizeLimitsPositive, fUseDirectionality,  fDirection_Filter, FilterThreshold);

	timer.Stop();
	fOutputProps.Vtx_Minimize_ComputeTime = timer.RealTime();

	fOutput.Vtx[0] = fRecoVtxPosFinal[0][0];
	fOutput.Vtx[1] = fRecoVtxPosFinal[0][1];
	fOutput.Vtx[2] = fRecoVtxPosFinal[0][2];
	fOutput.Vtx[3] = fRecoVtxPosFinal[0][3];
	fOutput.NLL = fRecoVtxPosFinal[0][4];
}

std::vector<VtxCandidate> LEAF::SearchVertex_Main(int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality)
{
// 1) Launch threads just like in SearchVertex_Main.
	std::vector<std::thread> lThreadList;
	std::vector<std::vector<double>> lOutputFinal;
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	int iCand_Step = (int)fPositionList.size()/fThread;
	if(iCand_Step < 1) iCand_Step = 1;

	for(int iStart=0; iStart<(int)fPositionList.size(); iStart+=iCand_Step)
	{
		std::thread tThrd(
			&LEAF::SearchVertex_thread,
			LEAF::GetME(),
			iStart, iCand_Step,
			nhits, tolerance,
			likelihood, lowerLimit, upperLimit, directionality
		);
		lThreadList.push_back(std::move(tThrd));
	}

	// std::cout << "time after vtx search threads in SearchVertex_Main : " << timer.RealTime() << std::endl;

	for(auto &th: lThreadList) if(th.joinable()) th.join();

	// 2) collect partial results from fThreadOutput.
	mtx.lock();
	lOutputFinal = fThreadOutput; // each row: [x, y, z, t, NLL]
	fThreadOutput.clear();
	mtx.unlock();

	// 3) Sort by NLL and keep 'tolerance' number of candidates
	std::sort(lOutputFinal.begin(), lOutputFinal.end(), SortOutputVector);
	while((int)lOutputFinal.size() > tolerance) lOutputFinal.pop_back();

	//4) Now build a vector of VtxCandidate.
	std::vector<VtxCandidate> finalList;
	finalList.reserve(lOutputFinal.size());

	for(auto &cand : lOutputFinal)
	{
		VtxCandidate cOut = CreateCandidate(cand, lowerLimit, upperLimit, false);
		cOut.NLL = cand[4]; // store the NLL from the coarse search.
		finalList.push_back(cOut);
	}

	return finalList;
}

// Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
// The step between each grid points is given by stepSize (in centimeters).
// We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
// The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower).
// The number of output vertices can be chosen using "tolerance".
void LEAF::SearchVertex_thread(int iStart, int iIte, int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality)
{
	std::vector<struct FitPosition> tVtxContainer;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	if (fPositionList.size() == 0)
		std::cout << " Position list not filled " << std::endl;
	if (VERBOSE >= 3)
		std::cout << "Number of hits for coarse grid search = " << nhits << std::endl;
	int iEnd = (iIte + iStart);
	if (iEnd > (int)fPositionList.size())
		iEnd = fPositionList.size();

	for (int iPos = iStart; iPos < iEnd; iPos++)
	{

		struct FitPosition hPos;
		hPos.Vtx = fPositionList[iPos];
		hPos.NLL = 0;

		if (likelihood) hPos.NLL = Likelihoods::Vertex_Time_NLL(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);
		else hPos.NLL = Likelihoods::Vertex_Score(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);

		tVtxContainer.push_back(hPos);
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), Likelihoods::SortingNLL());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();

	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	for (unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++)
	{
		std::vector<double> vPos(5, 0.);
		struct FitPosition hPos = tVtxContainer[iPos];

		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;
		fThreadOutput.push_back(vPos);
	}
	mtx.unlock();
}

std::vector<std::vector<double>> LEAF::MinimizeVertex_Main(std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality, std::vector<double>* fDirection_Filter, float FilterThreshold)
{
	std::vector<std::thread> lThreadList;
	std::vector<std::vector<double>> lOutputFinal;

	int iCand_Step = nCandidates / fThread;
	if (iCand_Step < 1) iCand_Step = 1;

	for (int iStart = 0; iStart < nCandidates; iStart += iCand_Step)
	{
		std::thread tThrd(&LEAF::MinimizeVertex_thread, LEAF::GetME(), iStart, iCand_Step,
						  initialVertex, limits, stepSize, nhits, nCandidates,
						  tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality, fDirection_Filter, FilterThreshold);

		lThreadList.push_back(std::move(tThrd));
	}

	for (unsigned int iIdx = 0; iIdx < lThreadList.size(); iIdx++)
	{
		// Wait thread
		if (lThreadList[iIdx].joinable())
			lThreadList[iIdx].join();
	}

	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	lOutputFinal = fThreadOutput;
	fThreadOutput.clear();
	mtx.unlock();

	std::sort(lOutputFinal.begin(), lOutputFinal.end(), SortOutputVector);

	while ((int)lOutputFinal.size() > tolerance) lOutputFinal.pop_back();

	return lOutputFinal;
}

// here is the function where minimization gonna take place. The loop is inside this function
void LEAF::MinimizeVertex_thread(
	int iStart, int iIte,
	std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits,
	int nCandidates, int tolerance, int verbose, bool /*likelihood*/, bool /*average*/, double lowerLimit, double upperLimit, int directionality, std::vector<double>* fDirection_Filter, float FilterThreshold)
{
	int iEnd = (iIte + iStart);
	if (iEnd > nCandidates) iEnd = nCandidates;

	std::vector<struct FitPosition> tVtxContainer;

	mtx.lock();
	TFitter *minimizer = new TFitter(14); // 4=nb de params?
	mtx.unlock();
	TMinuit *minuit = minimizer->GetMinuit();

	double arglist[20];
	int err = 0;
	// arglist[0]=0;
	double p1 = verbose - 1;
	minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1); // quiet mode
	minuit->SetErrorDef(1);
	// minuit->SetPrintLevel(verbose-1); // quiet mode
	if (VERBOSE < 2) minuit->mnexcm("SET NOWarnings", 0, 0, err);
	arglist[0] = 2;
	minuit->mnexcm("SET STR", arglist, 1, err); // set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
	minimizer->SetFCN(MinuitLikelihood);		// here is function to minimize

	double Mig[2] = {1e6, 1e0}; // maxcalls and tolerance

	minimizer->SetParameter(0, "vertex0", 0, stepSize, -limits[0], limits[0]);
	minimizer->SetParameter(1, "vertex1", 0, stepSize, -limits[1], limits[1]);
	minimizer->SetParameter(2, "vertex2", 0, stepSize, -limits[2], limits[2]);
	minimizer->SetParameter(3, "vertex3", 0, stepSize / fLightSpeed, -limits[3] / fLightSpeed, limits[3] / fLightSpeed);

	minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
	minimizer->SetParameter(5, "nhits", nhits, nhits, nhits - 1, nhits + 1);
	minimizer->SetParameter(6, "lowerLimit", lowerLimit, 5e-1, -7, -2);
	minimizer->SetParameter(7, "upperLimit", upperLimit, 5e-1, 2, 10);
	minimizer->SetParameter(8, "expoSigma", 100, 5, 0, 1e3);
	minimizer->SetParameter(9, "directionality", directionality, 1, 0, 2);
	minimizer->SetParameter(10, "dir_x", 0, 0, -1, 1);
	minimizer->SetParameter(11, "dir_y", 0, 0, -1, 1);
	minimizer->SetParameter(12, "dir_z", 0, 0, -1, 1);
	minimizer->SetParameter(13, "dir_threshold", 0, 0, 0, 0);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);

	for (int icand = iStart; icand < iEnd; icand++)
	{
		minimizer->SetParameter(0, "vertex0", initialVertex[icand][0], stepSize, initialVertex[icand][0] - limits[0], initialVertex[icand][0] + limits[0]);
		minimizer->SetParameter(1, "vertex1", initialVertex[icand][1], stepSize, initialVertex[icand][1] - limits[1], initialVertex[icand][1] + limits[1]);
		minimizer->SetParameter(2, "vertex2", initialVertex[icand][2], stepSize, initialVertex[icand][2] - limits[2], initialVertex[icand][2] + limits[2]);
		minimizer->SetParameter(3, "vertex3", initialVertex[icand][3], stepSize / fLightSpeed, initialVertex[icand][3] - limits[3] / fLightSpeed, initialVertex[icand][3] + limits[3] / fLightSpeed);

		minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
		minimizer->SetParameter(5, "nhits", nhits, nhits, nhits - 1, nhits + 1);
		minimizer->SetParameter(6, "lowerLimit", lowerLimit, 5e-1, -7, -2);
		minimizer->SetParameter(7, "upperLimit", upperLimit, 5e-1, 2, 10);
		minimizer->SetParameter(8, "expoSigma", 100, 5, 0, 1e3);
		minimizer->SetParameter(9, "directionality", directionality, 1, 0, 2);

		if(fDirection_Filter != nullptr && fDirection_Filter->size() > 2)
		{
			minimizer->SetParameter(10, "dir_x", (*fDirection_Filter)[0], 5e-1, (*fDirection_Filter)[0] - 1, (*fDirection_Filter)[0] + 1);
			minimizer->SetParameter(11, "dir_y", (*fDirection_Filter)[1], 5e-1, (*fDirection_Filter)[1] - 1, (*fDirection_Filter)[1] + 1);
			minimizer->SetParameter(12, "dir_z", (*fDirection_Filter)[2], 5e-1, (*fDirection_Filter)[2] - 1, (*fDirection_Filter)[2] + 1);
		}
		minimizer->SetParameter(13, "dir_threshold", FilterThreshold, 5e-1, FilterThreshold - 1, FilterThreshold + 1);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);
		minimizer->FixParameter(10);
		minimizer->FixParameter(11);
		minimizer->FixParameter(12);
		minimizer->FixParameter(13);

		minimizer->ExecuteCommand("MIGRAD", Mig, 2);

		std::vector<double> tCandRecoVtxPos(4, 0.);

		tCandRecoVtxPos[0] = minimizer->GetParameter(0);
		tCandRecoVtxPos[1] = minimizer->GetParameter(1);
		tCandRecoVtxPos[2] = minimizer->GetParameter(2);
		tCandRecoVtxPos[3] = minimizer->GetParameter(3);

		struct FitPosition hPos;
		hPos.Vtx = tCandRecoVtxPos;
		hPos.NLL = 0;

		// Get minuit
		minuit = minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);

		hPos.NLL = fmin;

		tVtxContainer.push_back(hPos);
		if (VERBOSE >= 2)
			std::cout << "Initial vertex position = " << initialVertex[icand][0] << "," << initialVertex[icand][1] << "," << initialVertex[icand][2] << "," << initialVertex[icand][3] << ", Candidate vertex position = (" << tCandRecoVtxPos[0] << "," << tCandRecoVtxPos[1] << "," << tCandRecoVtxPos[2] << "," << tCandRecoVtxPos[3] << ", NLL = " << fmin << std::endl;
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), Likelihoods::SortingNLL());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	for (unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++)
	{
		std::vector<double> vPos(5, 0.);
		struct FitPosition hPos = tVtxContainer[iPos];

		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;

		fThreadOutput.push_back(vPos);
	}
	delete minimizer;
	mtx.unlock();
}


/*****************************************************************************************************/
/* ENERGY FIT */
/*****************************************************************************************************/

void LEAF::FitEnergy()
{
	double totalPE = 0;  
	double totalPECorr = 0;

	TStopwatch timer;
	timer.Reset();
	timer.Start();
	
	// Taha's Polynomial Correction
	// double a = -0.000440808;
	// double b = 1.71313;
	// double c = -2.01675;
	// double d = 1.81755;
	
	for (unsigned int ihit = 0; ihit < fHitCollection->Size(); ihit++)
	{
		Hit lHit = fHitCollection->At(ihit);
		
		// int iPMT = lHit.PMT;
		// PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		// double PMTOrientation[3];
		// double dirVector[3];
		// for(int j = 0; j<3; j++) PMTOrientation[j] = lPMTInfo.Orientation[j];
		// for(int j = 0; j<3; j++) dirVector[j] = -(fTrueVtxPos[j] - lPMTInfo.Position[j]);
		// GeoTools::Normalize(PMTOrientation);
		// GeoTools::Normalize(dirVector);
		// double cosPMThitAngle = abs(dirVector[0] * PMTOrientation[0] + dirVector[1] * PMTOrientation[1] + dirVector[2] * PMTOrientation[2]);
		// double Ftheta = a + b*cosPMThitAngle + c*(pow(cosPMThitAngle,2)) + d*(pow(cosPMThitAngle,3)); // Taha's Polynomial Correction
		
		totalPECorr += lHit.Q;
		totalPE += lHit.Q;
	}
	
	double energyEstimate = totalPECorr * 0.111247 - 12.970173; //* Affine correction (affine because of dark rate and systematic errors), convert from charge to energy
	fOutput.TotalCharge = totalPE;
	fOutput.Energy = energyEstimate;

	timer.Stop();
	fOutputProps.Energy_Fit_ComputeTime = timer.RealTime();
}


/*****************************************************************************************************/
/* DIRECTION AND VERTEX FIT */
/*****************************************************************************************************/

std::vector<JointFitCandidate> LEAF::SearchVertexAndDir()
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<JointFitCandidate> Candidates = std::vector<JointFitCandidate>();

	std::vector<VtxCandidate> tVtxCandidates = this->SearchVertex_Main(fHitCollection->Size(),fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
	
	for(int i=0; i<(int)tVtxCandidates.size(); i++)
	{
		std::vector<double> vtx = {tVtxCandidates[i].X, tVtxCandidates[i].Y, tVtxCandidates[i].Z, tVtxCandidates[i].T};
		std::vector<DirectionCandidate> tDirCandidates = this->SearchDirection(vtx, DirTolerance);
		for(int j=0; j<(int)tDirCandidates.size(); j++)
		{
			JointFitCandidate candidate;
			candidate.VtxPart = tVtxCandidates[i];
			candidate.DirPart = tDirCandidates[j];
			candidate.NLL = candidate.DirPart.DNLL * candidate.VtxPart.NLL;
			Candidates.push_back(candidate);
		}
	}

	std::sort(Candidates.begin(), Candidates.end(), Likelihoods::CompareByNLL<JointFitCandidate>{});

	while ((int)Candidates.size() > fSearchVtxTolerance) Candidates.pop_back();

	std::vector<DirectionCandidate> dirCandidates = std::vector<DirectionCandidate>();
	for(int i=0; i<(int)Candidates.size(); i++) dirCandidates.push_back(Candidates[i].DirPart);
	// fOutput.stepOneContainsTrueDir = (int)ContainsTrueDir(&dirCandidates, fTrueDir);

	timer.Stop();
	fOutputProps.Vtx_Search_ComputeTime = timer.RealTime();

	return Candidates;
}

JointFitCandidate LEAF::MinimizeVertexAndDir(std::vector<JointFitCandidate> initialCandidates)
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	double stepSize = 10;
	int nhits = fHitCollection->Size();
	double *limits = new double[4];
	for (int i = 0; i < 4; i++) limits[i] = 2 * fSearchVtxStep;

	std::vector<struct JointFitCandidate> tContainer;

	TFitter *minimizer = new TFitter(12);
	TMinuit *minuit = minimizer->GetMinuit();

	double arglist[20];
	int err = 0;
	double p1 = VERBOSE - 1;
	minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1); // quiet mode
	minuit->SetErrorDef(1);
	if (VERBOSE < 2) minuit->mnexcm("SET NOWarnings", 0, 0, err);
	arglist[0] = 2;
	minuit->mnexcm("SET STR", arglist, 1, err); // set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
	minimizer->SetFCN(MinuitJointNLL);		// here is function to minimize

	double Mig[2] = {1e6, 1e0}; // maxcalls and tolerance

	minimizer->SetParameter(0, "vertex0", 0, stepSize, -limits[0], limits[0]);
	minimizer->SetParameter(1, "vertex1", 0, stepSize, -limits[1], limits[1]);
	minimizer->SetParameter(2, "vertex2", 0, stepSize, -limits[2], limits[2]);
	minimizer->SetParameter(3, "vertex3", 0, stepSize / fLightSpeed, -limits[3] / fLightSpeed, limits[3] / fLightSpeed);
	minimizer->SetParameter(4, "theta", 0, 2*M_PI/180, 0.0, TMath::Pi());
	minimizer->SetParameter(5, "phi", 0, 2*M_PI/180, -TMath::Pi(), TMath::Pi());

	minimizer->SetParameter(6, "PMTConfiguration", 0, 0, 0, 0); // Useless now
	minimizer->SetParameter(7, "nhits", nhits, nhits, nhits - 1, nhits + 1);
	minimizer->SetParameter(8, "lowerLimit", fSTimePDFLimitsQueueNegative, 5e-1, -7, -2);
	minimizer->SetParameter(9, "upperLimit", fSTimePDFLimitsQueuePositive, 5e-1, 2, 10);
	minimizer->SetParameter(10, "expoSigma", 100, 5, 0, 1e3);
	minimizer->SetParameter(11, "thread", 0, 0, 0, 0);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);
	minimizer->FixParameter(10);
	minimizer->FixParameter(11);

	for (int icand = 0; icand < (int)initialCandidates.size(); icand++)
	{
		minimizer->SetParameter(0, "vertex0", initialCandidates[icand].VtxPart.X, stepSize, initialCandidates[icand].VtxPart.X - limits[0], initialCandidates[icand].VtxPart.X + limits[0]);
		minimizer->SetParameter(1, "vertex1", initialCandidates[icand].VtxPart.Y, stepSize, initialCandidates[icand].VtxPart.Y - limits[1], initialCandidates[icand].VtxPart.Y + limits[1]);
		minimizer->SetParameter(2, "vertex2", initialCandidates[icand].VtxPart.Z, stepSize, initialCandidates[icand].VtxPart.Z - limits[2], initialCandidates[icand].VtxPart.Z + limits[2]);
		minimizer->SetParameter(3, "vertex3", initialCandidates[icand].VtxPart.T, stepSize / fLightSpeed, initialCandidates[icand].VtxPart.T - limits[3] / fLightSpeed, initialCandidates[icand].VtxPart.T + limits[3] / fLightSpeed);
		minimizer->SetParameter(4, "theta", initialCandidates[icand].DirPart.theta, 2*M_PI/180, 0.0, TMath::Pi());
		minimizer->SetParameter(5, "phi", initialCandidates[icand].DirPart.phi, 2*M_PI/180, -TMath::Pi(), TMath::Pi());

		minimizer->ExecuteCommand("MIGRAD", Mig, 2);

		JointFitCandidate tCandReco;

		tCandReco.VtxPart.X = minimizer->GetParameter(0);
		tCandReco.VtxPart.Y = minimizer->GetParameter(1);
		tCandReco.VtxPart.Z = minimizer->GetParameter(2);
		tCandReco.VtxPart.T = minimizer->GetParameter(3);
		tCandReco.DirPart.theta = minimizer->GetParameter(4);
		tCandReco.DirPart.phi = minimizer->GetParameter(5);

		// Get minuit
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);
		tCandReco.NLL = fmin;

		tContainer.push_back(tCandReco);
	}

	std::sort(tContainer.begin(), tContainer.end(), Likelihoods::CompareByNLL<JointFitCandidate>{});

	delete minimizer;

	timer.Stop();
	fOutputProps.Vtx_Minimize_ComputeTime = timer.RealTime();

	return tContainer[0];
}


