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

/*****************************************************************************************************/

// void MinimizeVertex_CallThread(int iStart, int iIte, std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality)
// {
// 	LEAF::GetME()->MinimizeVertex_thread(iStart, iIte, initialVertex, limits, stepSize, nhits, nCandidates, tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality);
// }

// void SearchVertex_CallThread(int iStart, int iIte, int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality)
// {
// 	LEAF::GetME()->SearchVertex_thread(iStart, iIte, nhits, tolerance, likelihood, lowerLimit, upperLimit, directionality);
// }
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

/*****************************************************************************************************/
void LEAF::Initialize(const Geometry *lGeometry)
{
	InitConfig(lGeometry);
	InitInputs(lGeometry);
	LoadSplines();
}

//Creates the candidates list after Coarse grid search, refines the fit, then outputs the best candidate
std::vector<double> LEAF::FitDirection(const std::vector<double>& fixedVertexPosition, FitterOutput &fOutput, int nhits) 
{
    // double theta_step = 5 * TMath::DegToRad(); 
    // double phi_step = 5 * TMath::DegToRad(); 
	double *tLimits = new double[4];  
	int verbose = VERBOSE;
	for (int i = 0; i < 4; i++) tLimits[i] = 2 * fSearchVtxStep;

    std::vector<DirectionCandidate> candidates = SearchDirection(fixedVertexPosition, DirTolerance);

	fOutput.stepOneContainsTrueDir = (int)ContainsTrueDir(&candidates, fTrueDir);

    if (verbose >= 2) 
	{
        for (unsigned int j = 0; j < candidates.size(); j++) std::cout << "Candidate theta = " << candidates[j].theta * TMath::RadToDeg() << " deg, phi = " << candidates[j].phi * TMath::RadToDeg() << " deg, NLL = " << candidates[j].DNLL << std::endl;
    }

	DirectionCandidate finalCandidate = MinimizeDirection(fixedVertexPosition, candidates, nhits, tLimits, verbose);

	fOutput.Dir = PolarToCartesianNorm({finalCandidate.theta, finalCandidate.phi});
	fOutput.DNLL = finalCandidate.DNLL;

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

		// Store fixed vertex position and number of hits in member variables
		// fFixedVertexPosition = fixedVertexPosition;
		// fNHits = nhits;

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

// Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
// The step between each grid points is given by stepSize (in centimeters).
// We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
// The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower).
// The number of output vertices can be chosen using "tolerance".
std::vector<std::vector<double>> LEAF::SearchVertex(int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality)
{
	// 2. How to set the search?
	// double stepSize = 0.5;//in m

	std::vector<double> tBestReconstructedVertexPosition;

	// double ** reconstructedVertexPosition = new double[4];

	std::vector<struct FitPosition> tVtxContainer;

	clock_t timeStart = clock();

	if (fPositionList.size() == 0)
		std::cout << " Position list not filled " << std::endl;
	if (VERBOSE >= 3)
		std::cout << "Number of hits for coarse grid search = " << nhits << std::endl;
	for (unsigned int iPos = 0; iPos < fPositionList.size(); iPos++)
	{
		struct FitPosition hPos;
		hPos.Vtx = fPositionList[iPos];
		hPos.NLL = 0;

		if (likelihood)
			hPos.NLL = Likelihoods::Vertex_Time_NLL(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);
		else
			hPos.NLL = Likelihoods::Vertex_Score(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);

		hPos.NLL += fRand->Uniform(0, 1e-5); // Here only not to have exactly the same value for 2NLL, and be removed from map

		tVtxContainer.push_back(hPos);

		if (VERBOSE >= 3)
			std::cout << "Test Pos: " << hPos.Vtx[3] << " " << hPos.Vtx[0] << " " << hPos.Vtx[1] << " " << hPos.Vtx[2] << " NLL = " << hPos.NLL << " " << likelihood << std::endl;

#ifdef VERBOSE_VTX
		// if(verbose && distanceToTrue < 2) std::cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<std::endl;
#endif
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), Likelihoods::SortingNLL());
	while ((int)tVtxContainer.size() > tolerance)
	{
#ifdef VERBOSE_VTX
		std::cout << " Remove end: " << tVtxContainer.back().NLL << " (first " << tVtxContainer[0].NLL << ") " << tVtxContainer.size() << std::endl;
#endif
		tVtxContainer.pop_back();
	}

	clock_t timeEndLoop = clock();

	std::vector<std::vector<double>> tReconstructedVertexPosition;
	for (int iPos = 0; iPos < tolerance; iPos++)
	{
		std::vector<double> vPos(5, 0.);
		struct FitPosition hPos = tVtxContainer[iPos];

		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;

		tReconstructedVertexPosition.push_back(vPos);

#ifdef VERBOSE_VTX
		std::cout << "NLL = " << hPos.NLL << ", vertex = " << vPos[0] << ", " << vPos[1] << ", " << vPos[2] << ", " << vPos[3] << std::endl;
#endif
	}

	clock_t timeEnd = clock();

	std::cout << "SearchVertex: Time it took for the loop = " << (timeEndLoop - timeStart) / 1e6 << ", time total = " << (timeEnd - timeStart) / 1e6 << " " << fPositionList.size() << std::endl;

	return tReconstructedVertexPosition;
}

std::vector<VtxCandidate> LEAF::SearchVertex_Main(int nhits, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality)
{
// 1) Launch threads just like in SearchVertex_Main.
	std::vector<std::thread> lThreadList;
	std::vector<std::vector<double>> lOutputFinal;

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
		// cand = [ x, y, z, t, NLL ]
		double vx = cand[0];
		double vy = cand[1];
		double vz = cand[2];
		double vt = cand[3];
		double nll= cand[4];

		// Build a 4-vector.
		std::vector<double> vtx(4);
		vtx[0] = vx;
		vtx[1] = vy;
		vtx[2] = vz;
		vtx[3] = vt;

		// 5) Per-candidate SNR.
		VtxCandidate cOut = ComputeCandidateSNR(vtx, lowerLimit, upperLimit);
		cOut.NLL = nll; // store the NLL from the coarse search.
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
	// 2. How to set the search?
	// double stepSize = 0.5;//in m

	std::vector<struct FitPosition> tVtxContainer;

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
		// std::cout << " Thread2 " << iPos << " " << iPos-iStart << " " << iEnd <<  std::endl;

		if (likelihood) hPos.NLL = Likelihoods::Vertex_Time_NLL(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);
		else hPos.NLL = Likelihoods::Vertex_Score(fHitCollection, hPos.Vtx, nhits, lowerLimit, upperLimit, false, false, directionality);

		tVtxContainer.push_back(hPos);
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), Likelihoods::SortingNLL());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();

	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	// std::vector< std::vector<double> > tReconstructedVertexPosition;
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

// Fine grid search of the vertex. It will search in a box this time, and not a cylinder as searchVertex is doing. It is useful when we have one or several first candidate vertices. Therefore, this code is generally run after serachVertex.
// The box is centered around the initial vertex provided.
// The box size is 2 x limits in each direction, centered on the initial vertex provided.
// nCandidate provide the size of the list of inital vertices, while tolerance provide the size of the output list.
// Other informations are the same as searchVertex.
std::vector<std::vector<double>> LEAF::SearchVertexFine(std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality)
{
	if (verbose) std::cout << "Start fine search, use directionality? " << directionality << std::endl;
	// 2. How to set the search?
	// double stepSize = 0.5;//in m
	// double * reconstructedVertexPosition = new double[4];
	std::vector<double> tBestReconstructedVertexPosition;
	double minNLL = 1e15;
	clock_t timeStart = clock();

	std::map<double, std::vector<double>> vertexContainer;

	for (int icand = 0; icand < nCandidates; icand++)
	{
		if (verbose)
			std::cout << "Candidate #" << icand << std::endl;
		for (double time = initialVertex[icand][VTX_T] - limits[0] / fLightSpeed; time <= initialVertex[icand][VTX_T] + limits[0] / fLightSpeed; time += (stepSize / fLightSpeed))
		{
			if (VERBOSE >= 2)
				std::cout << "Time = " << time << "ns" << std::endl;
			for (double x = initialVertex[icand][VTX_X] - limits[1]; x <= initialVertex[icand][VTX_X] + limits[1]; x += stepSize)
			{
				if (VERBOSE >= 2)
					std::cout << "x = " << x << "m" << std::endl;
				for (double y = initialVertex[icand][VTX_Y] - limits[2]; y <= initialVertex[icand][VTX_Y] + limits[2]; y += stepSize)
				{
					if (VERBOSE >= 2)
						std::cout << "y = " << y << "m" << std::endl;
					for (double z = initialVertex[icand][VTX_Z] - limits[3]; z <= initialVertex[icand][VTX_Z] + limits[3]; z += stepSize)
					{
						if (VERBOSE >= 2)
							std::cout << "z = " << z << "m" << std::endl;

						std::vector<double> tVtxPosition(4, 0.); // In centimeters
						tVtxPosition[0] = x;					 // radius*TMath::Cos(angle);
						tVtxPosition[1] = y;					 // radius*TMath::Sin(angle);
						tVtxPosition[2] = z;					 // height;
						tVtxPosition[3] = time;

						double NLL = 0;
						if (average)
						{
							// double vNLL[fAveraging];
							for (int irand = 0; irand < fAveraging; irand++)
							{
								// Determine the position around the vertex
								double radius = 1.5 * fRand->Uniform(0, stepSize);
								double phi = fRand->Uniform(0, 2 * TMath::Pi());
								double theta = fRand->Uniform(0, TMath::Pi());
								double timeRand = 1.5 * fRand->Uniform(-stepSize / fLightSpeed, stepSize / fLightSpeed);
								std::vector<double> randVertexPosition(4, 0.); // In centimeters
								randVertexPosition[0] = tVtxPosition[0] + radius * TMath::Sin(theta) * TMath::Cos(phi);
								randVertexPosition[1] = tVtxPosition[1] + radius * TMath::Sin(theta) * TMath::Sin(phi);
								randVertexPosition[2] = tVtxPosition[2] + radius * TMath::Cos(theta);
								randVertexPosition[3] = tVtxPosition[3] + timeRand; // rand->Uniform(-stepSize/fLightSpeed/2,stepSize/fLightSpeed/2);
								// Evaluate its NLL
								double nll = Likelihoods::FindNLL(fHitCollection, randVertexPosition, nhits, likelihood, verbose, lowerLimit, upperLimit, false, false, directionality);
								NLL += nll;
								// vNLL[irand]=nll;
							}

							if (fAveraging != 0)
								NLL /= fAveraging;
						}
						else
							NLL = Likelihoods::FindNLL(fHitCollection, tVtxPosition, nhits, likelihood, verbose, lowerLimit, upperLimit, false, false, directionality);

						std::vector<double> candidateReconstructedVertexPosition(4, 0.);
						for (int i = 0; i < 4; i++)
							candidateReconstructedVertexPosition[i] = tVtxPosition[i];

						vertexContainer[NLL] = candidateReconstructedVertexPosition;
						if ((int)vertexContainer.size() > tolerance)
						{
							vertexContainer.erase(--vertexContainer.rbegin().base());
						}

						if (NLL < minNLL)
						{
							minNLL = NLL;
							// std::cout << "minNLL = "<<minNLL << std::endl;
							tBestReconstructedVertexPosition = tVtxPosition;
						}
					}
				}
			}
		}
	}
	clock_t timeEndLoop = clock();
	std::vector<std::vector<double>> tReconstructedVertexPosition;
	for (int count = 0; count < tolerance; count++)
	{
		std::vector<double> vTmp(5, 0.);
		tReconstructedVertexPosition.push_back(vTmp);
	}
	int counter = 0;
	for (std::map<double, std::vector<double>>::iterator it = vertexContainer.begin(); it != vertexContainer.end(); it++)
	{
		// for(std::map< double, std::vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
		if (counter < tolerance)
		{
			// tReconstructedVertexPosition[counter] = new double[4];
			for (int i = 0; i < 4; i++)
				tReconstructedVertexPosition[counter][i] = (it->second)[i];
			tReconstructedVertexPosition[counter][4] = it->first;
		}
		counter++;
	}

	clock_t timeEnd = clock();
	if (verbose) std::cout << "SearchVertexFine: Time it took for the loop = " << (timeEndLoop - timeStart) / 1e6 << ", time total = " << (timeEnd - timeStart) / 1e6 << std::endl << std::endl;

	return tReconstructedVertexPosition;
	// return tBestReconstructedVertexPosition;
}

// here is the function where minimization gonna take place. The loop is inside this function
std::vector<std::vector<double>> LEAF::MinimizeVertex(std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits, int nCandidates, int tolerance, int verbose, bool /*likelihood*/, bool /*average*/, double lowerLimit, double upperLimit, int directionality)
{
	if (VERBOSE >= 2) std::cout << "Minimizer" << std::endl;
	// 2. How to set the search?
	// double stepSize = 0.5;//in m
	// double * reconstructedVertexPosition = new double[4];
	std::vector<double> tBestReconstructedVertexPosition;
	clock_t timeStart = clock();

	std::vector<struct FitPosition> tVtxContainer;

	if (VERBOSE >= 2)
		std::cout << "Minimizer" << std::endl;
	TFitter *minimizer = new TFitter(8); // 4=nb de params?
	TMinuit *minuit = minimizer->GetMinuit();

	double arglist[20];
	int err = 0;
	// arglist[0]=0;
	double p1 = verbose - 1;
	minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1); // quiet mode
	minuit->SetErrorDef(1);
	// minuit->SetPrintLevel(-1); // quiet mode
	if (VERBOSE < 2) minuit->mnexcm("SET NOWarnings", 0, 0, err);
	arglist[0] = 2;
	minuit->mnexcm("SET STR", arglist, 1, err); // set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
	minimizer->SetFCN(MinuitLikelihood);		// here is function to minimize

	double Mig[2] = {1e6, 1e0}; // maxcalls and tolerance
	// double Imp2[1]={10000};

	minimizer->SetParameter(0, "vertex0", 0, stepSize, -limits[0], limits[0]);
	minimizer->SetParameter(1, "vertex1", 0, stepSize, -limits[1], limits[1]);
	minimizer->SetParameter(2, "vertex2", 0, stepSize, -limits[2], limits[2]);
	minimizer->SetParameter(3, "vertex3", 0, stepSize / fLightSpeed, -limits[3] / fLightSpeed, limits[3] / fLightSpeed);

	minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
	minimizer->SetParameter(5, "nhits", nhits, nhits, nhits - 1, nhits + 1);
	minimizer->SetParameter(6, "lowerLimit", lowerLimit, 5e-1, -7, -2);
	minimizer->SetParameter(7, "upperLimit", upperLimit, 5e-1, 2, 10);
	// minimizer->SetParameter(8,"scaling factor",1,1,0,1e3);
	// minimizer->SetParameter(8,"signal pe",1,1,0,1e3);
	minimizer->SetParameter(8, "expoSigma", 100, 5, 0, 1e3);
	minimizer->SetParameter(9, "directionality", directionality, 1, 0, 2);
	minimizer->SetParameter(10, "thread", 0, 0, 0, 0);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);
	minimizer->FixParameter(10);

	for (int icand = 0; icand < nCandidates; icand++)
	{
#ifdef VERBOSE_MIN
		if (verbose >= 1)
		{
			std::cout << "Candidate #" << icand << std::endl;
			std::cout << "Initial (x,y,z,t): " << initialVertex[icand][0]
					  << ", " << initialVertex[icand][1]
					  << ", " << initialVertex[icand][2]
					  << ", " << initialVertex[icand][3] << std::endl;
		}
#endif

		minimizer->SetParameter(0, "vertex0", initialVertex[icand][0], stepSize, initialVertex[icand][0] - limits[0], initialVertex[icand][0] + limits[0]);
		minimizer->SetParameter(1, "vertex1", initialVertex[icand][1], stepSize, initialVertex[icand][1] - limits[1], initialVertex[icand][1] + limits[1]);
		minimizer->SetParameter(2, "vertex2", initialVertex[icand][2], stepSize, initialVertex[icand][2] - limits[2], initialVertex[icand][2] + limits[2]);
		minimizer->SetParameter(3, "vertex3", initialVertex[icand][3], stepSize / fLightSpeed, initialVertex[icand][3] - limits[3] / fLightSpeed, initialVertex[icand][3] + limits[3] / fLightSpeed);

		minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
		minimizer->SetParameter(5, "nhits", nhits, nhits, nhits - 1, nhits + 1);
		minimizer->SetParameter(6, "lowerLimit", lowerLimit, 5e-1, -7, -2);
		minimizer->SetParameter(7, "upperLimit", upperLimit, 5e-1, 2, 10);
		minimizer->SetParameter(8, "expoSigma", 100, 5, 0, 1e3);
		// minimizer->SetParameter(8,"scaling factor",1,1,0,1e3);
		// minimizer->SetParameter(8,"signal pe",1,1,0,1e3);
		minimizer->SetParameter(9, "directionality", directionality, 1, 0, 2);
		minimizer->SetParameter(10, "thread", 0, 0, 0, 0);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);
		minimizer->FixParameter(10);

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

#ifdef VERBOSE_MIN
		if (verbose >= 1)
			std::cout << "NLL=" << hPos.NLL << std::endl;
#endif
		tVtxContainer.push_back(hPos);

		if (VERBOSE >= 2)
			std::cout << "Initial vertex position = " << initialVertex[icand][0] << "," << initialVertex[icand][1] << "," << initialVertex[icand][2] << "," << initialVertex[icand][3] << ", Candidate vertex position = (" << tCandRecoVtxPos[0] << "," << tCandRecoVtxPos[1] << "," << tCandRecoVtxPos[2] << "," << tCandRecoVtxPos[3] << ", NLL = " << fmin << std::endl;
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), Likelihoods::SortingNLL());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();

	clock_t timeEndLoop = clock();

	std::vector<std::vector<double>> tReconstructedVertexPosition;
	for (int iPos = 0; iPos < tolerance; iPos++)
	{
		std::vector<double> vPos(5, 0.);
		struct FitPosition hPos = tVtxContainer[iPos];

		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;

		tReconstructedVertexPosition.push_back(vPos);
	}

	delete minimizer;

	tBestReconstructedVertexPosition = tReconstructedVertexPosition[0];

	clock_t timeEnd = clock();

	if (verbose) std::cout << "MinimizeVertex: Time it took for the loop = " << (timeEndLoop - timeStart) / 1e6 << ", time total = " << (timeEnd - timeStart) / 1e6 << std::endl << std::endl;

	return tReconstructedVertexPosition;
	// return tBestReconstructedVertexPosition;
}

std::vector<std::vector<double>> LEAF::MinimizeVertex_Main(std::vector<std::vector<double>> initialVertex, double *limits, double stepSize, int nhits, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality)
{
	std::vector<std::thread> lThreadList;
	std::vector<std::vector<double>> lOutputFinal;

	int iCand_Step = nCandidates / fThread;
	if (iCand_Step < 1) iCand_Step = 1;

	for (int iStart = 0; iStart < nCandidates; iStart += iCand_Step)
	{
		// std::thread tThrd(MinimizeVertex_CallThread, iStart, iCand_Step,
		// 				  initialVertex, limits, stepSize, nhits, nCandidates,
		// 				  tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality);

		std::thread tThrd(&LEAF::MinimizeVertex_thread, LEAF::GetME(), iStart, iCand_Step,
						  initialVertex, limits, stepSize, nhits, nCandidates,
						  tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality);

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
	int nCandidates, int tolerance, int verbose, bool /*likelihood*/, bool /*average*/, double lowerLimit, double upperLimit, int directionality)
{
	int iEnd = (iIte + iStart);
	if (iEnd > nCandidates) iEnd = nCandidates;

	std::vector<struct FitPosition> tVtxContainer;

	mtx.lock();
	TFitter *minimizer = new TFitter(8); // 4=nb de params?
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
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);

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
	while ((int)tVtxContainer.size() > tolerance) 
	{
		// std::cout << " Remove end: " << tVtxContainer.back().NLL << " (first " << tVtxContainer[0].NLL << ") " << tVtxContainer.size() << std::endl;
		tVtxContainer.pop_back();
	}

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

void LEAF::FitVertex(FitterOutput &fOutput)
{
	double *tLimits = new double[4];
	int iHitsTotal = fHitCollection->Size();
	std::vector< std::vector<double> > fRecoVtxPosFinal;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<VtxCandidate> tRecoVtxPosCand = this->SearchVertex_Main(iHitsTotal,fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
	std::vector<std::vector<double>> tRecoVtxPos;
	std::vector<double> SNRList; 

	// Extract only the vertex positions from the VtxCandidate structure
	for (const auto& cand : tRecoVtxPosCand) 
	{
		tRecoVtxPos.push_back({cand.X, cand.Y, cand.Z, cand.T, cand.NLL});
		SNRList.push_back(cand.SNR); // Store corresponding SNR
	}

	std::vector<double> trueVertexPosition(4);
	trueVertexPosition[0] = fTrueVtxPos[0]; // x
	trueVertexPosition[1] = fTrueVtxPos[1]; // y
	trueVertexPosition[2] = fTrueVtxPos[2]; // z
	trueVertexPosition[3] = fTrueVtxPos[3]; // t


	if(ContainsTrueVtx(&tRecoVtxPosCand, fTrueVtxPos))
		fOutput.stepOneContainsTrueVtx = 1;

	if (VERBOSE >= 2)
	{
		for (unsigned int a = 0; a < tRecoVtxPos.size(); a++)
			std::cout << "Candidate vertex [ " << a << " ] = " << tRecoVtxPos[a][0] << " , " << tRecoVtxPos[a][1] << " , " << tRecoVtxPos[a][2] << " , " << tRecoVtxPos[a][3] << ", NLL = " << tRecoVtxPos[a][4] << std::endl;
	}

	timer.Stop();
	fOutput.firstStepTime = timer.RealTime();
	std::cout << "SearchVertex took: " << timer.RealTime() << std::endl;

	double dStepSizeFinal = 10;
	int iToleranceFinal = 1;
	if (fStepByStep) //* Not used and not tested
	{
		// 2.b. Refined vertex search around candidate found previously
		//  -> Commented since ages, removed on 2020/11/20

		// 2.c. More refined vertex search around candidate found previously
		double dStepSizeFine2 = 100;
		int iToleranceFine2 = 20;
		for (int i = 0; i < 4; i++) tLimits[i] = fSearchVtxStep; // particleStart[i]*1e-2;//onversion to meter
		std::vector<std::vector<double>> tRecoVtxPosFine2 = this->SearchVertexFine(tRecoVtxPos, tLimits, dStepSizeFine2, iHitsTotal, fSearchVtxTolerance, iToleranceFine2, VERBOSE, true, false, fSTimePDFLimitsQueueNegative, fSTimePDFLimitsQueuePositive, fUseDirectionality);
		//////////////////////////////////////////

		// 2.d. Final refined vertex search around candidate found previously
		double dStepSizeFine3 = 20;
		int iToleranceFine3 = 1;
		for (int i = 0; i < 4; i++)
			tLimits[i] = dStepSizeFine2; // particleStart[i]*1e-2;//onversion to meter
		std::vector<std::vector<double>> tRecoVtxPosFine3 = this->SearchVertexFine(tRecoVtxPosFine2, tLimits, dStepSizeFine3, iHitsTotal, iToleranceFine2, iToleranceFine3, VERBOSE, true, false, fSTimePDFLimitsQueueNegative, fSTimePDFLimitsQueuePositive, fUseDirectionality);
		if (VERBOSE == 1)
		{
			for (int a = 0; a < iToleranceFine3; a++)
				std::cout << "Final candidate vertex time = " << tRecoVtxPosFine3[a][0] << ", x = " << tRecoVtxPosFine3[a][1] << ", y = " << tRecoVtxPosFine3[a][2] << ", z=" << tRecoVtxPosFine3[a][3] << std::endl;
		}
		//////////////////////////////////////////

		for (int i = 0; i < 4; i++)
			tLimits[i] = 100; // dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
		fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPosFine3, tLimits, dStepSizeFinal, iHitsTotal, iToleranceFine3, iToleranceFinal, VERBOSE, true, false, fMinimizeLimitsNegative, fMinimizeLimitsPositive, fUseDirectionality);
		if (VERBOSE >= 1)
		{
			for (int a = 0; a < iToleranceFinal; a++)
			{
				std::cout << "Final candidate vertex time = " << fRecoVtxPosFinal[a][0] << ", x = " << fRecoVtxPosFinal[a][1] << ", y = " << fRecoVtxPosFinal[a][2] << ", z=" << fRecoVtxPosFinal[a][3] << std::endl;
			}
		}
	}
	else
	{
		timer.Reset();
		timer.Start();
		
		for (int i = 0; i < 4; i++) tLimits[i] = 2 * fSearchVtxStep;

		// fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPos,tLimits,dStepSizeFinal,iHitsTotal,fSearchVtxTolerance,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,fUseDirectionality);
		fRecoVtxPosFinal = this->MinimizeVertex_Main(tRecoVtxPos, tLimits, dStepSizeFinal, iHitsTotal, fSearchVtxTolerance, iToleranceFinal, VERBOSE, true, false, fMinimizeLimitsNegative, fMinimizeLimitsPositive, fUseDirectionality);

		timer.Stop();
		fOutput.secondStepTime = timer.RealTime();
		std::cout << "Minimizer took: " << timer.RealTime() << std::endl;
	}

	fOutput.Vtx[0] = fRecoVtxPosFinal[0][0];
	fOutput.Vtx[1] = fRecoVtxPosFinal[0][1];
	fOutput.Vtx[2] = fRecoVtxPosFinal[0][2];
	fOutput.Vtx[3] = fRecoVtxPosFinal[0][3];
	fOutput.NLL = fRecoVtxPosFinal[0][4];
	// fOutput->NLLR = Likelihoods::GoodnessOfFit(fHitCollection, fRecoVtxPosFinal[0], iHitsTotal, fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality );
	fOutput.SNRList = SNRList;
}

void LEAF::FitEnergy(FitterOutput &fOutput)
{
	int iInTime = 0;
	double totalPE = 0;  
	double totalPECorr = 0;
	
	double a = -0.000440808;
	double b = 1.71313;
	double c = -2.01675;
	double d = 1.81755;
	
	for (unsigned int ihit = 0; ihit < fHitCollection->Size(); ihit++)
	{
		// Hit lHit = fHitInfo[ihit];
		Hit lHit = fHitCollection->At(ihit);
		
		int iPMT = lHit.PMT;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		
		double PMTOrientation[3];
		for(int j = 0; j<3; j++) PMTOrientation[j] = lPMTInfo.Orientation[j];
		double dirVector[3];
		dirVector[0] = -(fTrueVtxPos[0] - lPMTInfo.Position[0]);
		dirVector[1] = -(fTrueVtxPos[1] - lPMTInfo.Position[1]);
		dirVector[2] = -(fTrueVtxPos[2] - lPMTInfo.Position[2]);
		GeoTools::Normalize(PMTOrientation);
		GeoTools::Normalize(dirVector);
		double cosPMThitAngle = abs(dirVector[0] * PMTOrientation[0] + dirVector[1] * PMTOrientation[1] + dirVector[2] * PMTOrientation[2]);
		double Ftheta = a + b*cosPMThitAngle + c*(pow(cosPMThitAngle,2)) + d*(pow(cosPMThitAngle,3));
		
		totalPECorr += lHit.Q*Ftheta/(cosPMThitAngle/**sin(PMThitAngle)*/); 
		totalPE += lHit.Q;
		
		// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		// double distance = Astro_GetDistance(lPMTInfo.Position, fOutput.Vtx);
		// double tof = distance / fLightSpeed;
		// double residual = hitTime - tof - fOutput.Vtx[3];
		double residual = ComputeResidualTime(fOutput.Vtx, fOutput.Vtx[3], fHitCollection->At(ihit));
		if (residual > fSTimePDFLimitsQueueNegative && residual < fSTimePDFLimitsQueuePositive) iInTime++;
	}
	
	fOutput.InTime 	= iInTime;
	double energyEstimate = totalPECorr ;
	fOutput.TotalCharge = totalPE;
	fOutput.Energy = energyEstimate;

}

std::vector<JointFitCandidate> LEAF::SearchVertexAndDir(FitterOutput &fOutput)
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<JointFitCandidate> Candidates = std::vector<JointFitCandidate>();

	std::vector<VtxCandidate> tVtxCandidates = this->SearchVertex_Main(fHitCollection->Size(),fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
	
	fOutput.stepOneContainsTrueVtx = (int)ContainsTrueVtx(&tVtxCandidates, fTrueVtxPos);
	
	for(int i=0; i<(int)tVtxCandidates.size(); i++)
	{
		std::vector<double> vtx = {tVtxCandidates[i].X, tVtxCandidates[i].Y, tVtxCandidates[i].Z, tVtxCandidates[i].T};
		std::vector<DirectionCandidate> tDirCandidates = this->SearchDirection(vtx, DirTolerance);
		for(int j=0; j<(int)tDirCandidates.size(); j++)
		{
			JointFitCandidate candidate;
			candidate.VtxPart = tVtxCandidates[i];
			candidate.DirPart = tDirCandidates[j];
			// double DNLL = Likelihoods::Dir_NLL(fHitCollection ,candidate.Vtx, theta, phi, fHitCollection->Size());
			// double VNLL = Likelihoods::Vertex_Score(fHitCollection, candidate.Vtx, fHitCollection->Size(), fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality);
			candidate.NLL = candidate.DirPart.DNLL * candidate.VtxPart.NLL;
			Candidates.push_back(candidate);
		}
	}

	std::sort(Candidates.begin(), Candidates.end(), Likelihoods::CompareByNLL<JointFitCandidate>{});

	while ((int)Candidates.size() > fSearchVtxTolerance) Candidates.pop_back();

	std::vector<DirectionCandidate> dirCandidates = std::vector<DirectionCandidate>();
	for(int i=0; i<(int)Candidates.size(); i++) dirCandidates.push_back(Candidates[i].DirPart);
	fOutput.stepOneContainsTrueDir = (int)ContainsTrueDir(&dirCandidates, fTrueDir);

	timer.Stop();
	fOutput.firstStepTime = timer.RealTime();

	return Candidates;
}

JointFitCandidate LEAF::MinimizeVertexAndDir(std::vector<JointFitCandidate> initialCandidates, FitterOutput &fOutput)
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
	fOutput.secondStepTime = timer.RealTime();

	return tContainer[0];
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
	fOutput.SNRList = std::vector<double>(); 
	fOutput.Vtx = std::vector<double>(4);
	fOutput.Dir = std::vector<double>(3);
	fOutput.NLL		= 0.;
	fOutput.Energy = 0.;
	fOutput.TotalCharge = 0.;
	fOutput.InTime		= 0;
	fOutput.True_NLLDiff	= 0;
	fOutput.True_TimeDiff	= 0;
	fOutput.True_TistDiff	= 0;
	return fOutput;
}

struct FitterOutput LEAF::MakeSequentialFit(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT)
{
	LoadHitCollection(lHitCol, lTriggerTime, bMultiPMT);
	
	struct FitterOutput fOutput = NewOutuput();
	if(fHitCollection->Size() <= 0) return fOutput;
	
	//* Vertex
	FitVertex(fOutput);

	//* Goodness of fit
	fOutput.NLLR = Likelihoods::GoodnessOfFit(fHitCollection, fOutput.Vtx, fHitCollection->Size(), fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality); //* Goodness of Fit
	
	//* Direction
	FitDirection(fOutput.Vtx, fOutput, fHitCollection->Size());

	//* Energy
	FitEnergy(fOutput);

	return fOutput;
}

struct FitterOutput LEAF::MakeJointFit(const HitCollection<Hit> *lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT)
{
	LoadHitCollection(lHitCol, lTriggerTime, bMultiPMT);

	struct FitterOutput fOutput = NewOutuput();
	if(fHitCollection->Size() <= 0) return fOutput;

	std::vector<JointFitCandidate> Candidates = SearchVertexAndDir(fOutput);
	JointFitCandidate finalCandidate = MinimizeVertexAndDir(Candidates, fOutput);

	fOutput.Vtx = {finalCandidate.VtxPart.X, finalCandidate.VtxPart.Y, finalCandidate.VtxPart.Z, finalCandidate.VtxPart.T};
	fOutput.Dir = PolarToCartesianNorm({finalCandidate.DirPart.theta, finalCandidate.DirPart.phi});
	fOutput.NLL = finalCandidate.NLL;
	fOutput.DNLL = Likelihoods::Dir_NLL(fHitCollection, fOutput.Vtx, finalCandidate.DirPart.theta, finalCandidate.DirPart.phi , fHitCollection->Size());

	fOutput.NLLR = Likelihoods::GoodnessOfFit(fHitCollection, fOutput.Vtx, fHitCollection->Size(), fMinimizeLimitsNegative, fMinimizeLimitsPositive, true, false, fUseDirectionality); //* Goodness of Fit
	
	//* Energy
	FitEnergy(fOutput);

	return fOutput;
}
