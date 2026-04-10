/*****************************************************************************************************/
/**	LEAFSrcFitter.cc																				**/
/**	Date: March 27th 2026																			**/
/**	Desc: Implementation of the Fitter functions													**/
/*****************************************************************************************************/

#ifndef LEAF_HPP
#include "LEAF.hpp"
#endif

/*****************************************************************************************************/
/* DIRECTION FIT */
/*****************************************************************************************************/

//Creates the candidates list after Coarse grid search, refines the fit, then outputs the best candidate
ROOT::Math::XYZVector LEAF::FitDirection(const ROOT::Math::XYZTVector& vertex, bool searchPrior) {
	//* Init
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	ROOT::Math::XYZTVector limits(1.,1.,1.,1.); 
	limits *= 2.*LEAFConfig::fSearchVtxStep;
	
	int verbose = VERBOSE;

	//* Fit first step : find candidates for the direction
	std::vector<DirectionCandidate> candidates;
	if(searchPrior) {
		//* Only one candidate, chosen wisely by computing a mean hit direction for hits within a specified residual time window
		ROOT::Math::XYZVector dirPriorCart = this->FitDirectionQuick(vertex);
		ROOT::Math::Polar3DVector dirPriorPol(dirPriorCart);

		DirectionCandidate priorCandidate;
		priorCandidate.theta = dirPriorPol.Theta();
		priorCandidate.phi = dirPriorPol.Phi();
		candidates.push_back(priorCandidate);
	}
	else {
		candidates = this->SearchDirection(vertex, LEAFConfig::fDirTolerance); //* multiple candidates, spread uniformly
	}

	// std::cout << "nb direction candidates : " << candidates.size() << std::endl;

	fPerformances.dir_search_ct = timer.RealTime();
	TStopwatch timer2;
	timer2.Reset();
	timer2.Start();

	//? To be removed for real data where true dir is not known, no impact on the fit
	// fOutput.stepOneContainsTrueDir = (int)ContainsTrueDir(&candidates, fTrueDir); //* a quick check to see if the step one works (one of the candidate is near the true direction)

	if (verbose >= 2) {
		for (unsigned int j = 0; j < candidates.size(); j++) {
			std::cout 	<< "Candidate "
						<< "theta = " << candidates[j].theta * TMath::RadToDeg() << " deg, "
						<< "phi = "   << candidates[j].phi   * TMath::RadToDeg() << " deg, "
						<< "NLL = " << candidates[j].nll 
						<< std::endl;
		} 
	}

	DirectionCandidate finalCandidate = this->MinimizeDirection(vertex, candidates, limits, verbose);

	ROOT::Math::Polar3DVector dirPol(1, finalCandidate.theta, finalCandidate.phi);
	fOutput.dir = ROOT::Math::XYZVector(dirPol.Unit()); // Unit is probably not needed?
	fOutput.dir_nll = finalCandidate.nll;

	fPerformances.dir_minimize_ct = timer2.RealTime();

	return fOutput.dir;
}

ROOT::Math::XYZVector LEAF::FitDirectionQuick(const ROOT::Math::XYZTVector& vertex) {
	// Quick fit implementation

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());

	std::vector<ROOT::Math::XYZVector> scaledDistToPMTs;
	std::vector<double> residuals;

	double dir_x = 0.0;
	double dir_y = 0.0;
	double dir_z = 0.0;

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for(unsigned int i = 0; i < fHitsCollection[iPMTType]->size(); i++) {
			std::vector<double> toHit;
			
			const HKHit* lHit = fHitsCollection[iPMTType]->at(i);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			ROOT::Math::XYZVector distToPMT = (lPMTInfo->GetPositionInCm() - vertexPos);
			double residual = LEAFUtilities::GetResidual(distToPMT, lHit->GetTime(), vertex.T());

			if ( residual < LEAFConfig::fDirectionPDF_maxResidual && residual > LEAFConfig::fDirectionPDF_minResidual) {
				ROOT::Math::XYZVector scaledDist = (distToPMT.Unit() * lHit->GetCharge());
				dir_x += scaledDist.X();
				dir_y += scaledDist.Y();
				dir_z += scaledDist.Z();
			}
		}
	}

	fPerformances.dir_quick_search_ct = timer.RealTime();
	fOutput.dir_quick = ROOT::Math::XYZVector(dir_x, dir_y, dir_z).Unit();
	return ROOT::Math::XYZVector(dir_x, dir_y, dir_z).Unit();
}

std::vector<LEAF::DirectionCandidate> LEAF::SearchDirection(const ROOT::Math::XYZTVector& vertex, int tolerance)
{
	std::vector<DirectionCandidate> candidates;

	// Grid search over theta and phi
	for (double theta = 0.0; theta <= TMath::Pi(); theta += LEAFConfig::fThetaStep) {
		for (double phi = -TMath::Pi(); phi < TMath::Pi(); phi += LEAFConfig::fPhiStep) {
			// Store the candidate
			DirectionCandidate candidate;
			candidate.theta = theta;
			candidate.phi = phi;
			candidate.nll = this->Dir_NLL(vertex, theta, phi);
			candidates.push_back(candidate);
		}
	}


	std::sort(candidates.begin(), candidates.end(), LEAF::SortingNLL<DirectionCandidate>());

	// Select the best candidates based on tolerance
	while((int)candidates.size() > tolerance) candidates.pop_back();

	return candidates;
}

//MIGRAD Optimization
LEAF::DirectionCandidate LEAF::MinimizeDirection(const ROOT::Math::XYZTVector& vertex, std::vector<DirectionCandidate>& candidates, ROOT::Math::XYZTVector limits, int verbose) {
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
	minimizer->SetFCN(LEAFLikelihoods::MinuitDirNLL);

	for (int icand = 0; icand < (int)candidates.size(); icand++) {
		// Set initial parameters and bounds
		minimizer->SetParameter(0, "theta", candidates[icand].theta, 2*TMath::Pi()/180, 0.0, TMath::Pi());
		minimizer->SetParameter(1, "phi", candidates[icand].phi, 2*TMath::Pi()/180, -TMath::Pi(), TMath::Pi());
		minimizer->SetParameter(2, "nhits", 0, 0, -1, +1);
		minimizer->SetParameter(3, "vertex0", vertex.X(), 5e-1, -limits.X(), limits.X());
		minimizer->SetParameter(4, "vertex1", vertex.Y(), 5e-1, -limits.Y(), limits.Y());
		minimizer->SetParameter(5, "vertex2", vertex.Z(), 5e-1, -limits.Z(), limits.Z());
		minimizer->SetParameter(6, "vertex3", vertex.T(), 5e-1 / LEAFConfig::fLightSpeed, -limits.T() / LEAFConfig::fLightSpeed, limits.T() / LEAFConfig::fLightSpeed);
		
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
		candidates[icand].theta = minimizer->GetParameter(0); // theta
		candidates[icand].phi = minimizer->GetParameter(1); // phi

		minuit = minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);

		candidates[icand].nll = fmin;
	}

	std::sort(candidates.begin(), candidates.end(), LEAF::SortingNLL<DirectionCandidate>());
	// Optional: Print the refined parameters
	// if (verbose >= 2) std::cout << "After refinement: theta = " << optimizedParameters[0] * TMath::RadToDeg() << " deg, phi = " << optimizedParameters[1] * TMath::RadToDeg() << " deg" << std::endl;

	delete minimizer;
	return candidates[0]; // Returns [theta, phi]
}


/*****************************************************************************************************/
/* VERTEX FIT */
/*****************************************************************************************************/


void LEAF::FitVertex(float FilterThreshold) {
	std::vector<LEAFStructures::VtxCandidate> fRecoVtxPosFinal;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<LEAFStructures::VtxCandidate> tRecoVtxPos = this->SearchVertex_Main(LEAFConfig::fSearchVtxTolerance,false,LEAFConfig::fSTimePDFLimitsQueueNegative,LEAFConfig::fSTimePDFLimitsQueuePositive,false);
	
	if (VERBOSE >= 2) {
		for (unsigned int i=0; i<tRecoVtxPos.size(); i++) {
			std::cout << "Candidate vertex [ " << i << " ] = (" 
					  << tRecoVtxPos[i].vtx.X() << " , " 
					  << tRecoVtxPos[i].vtx.Y() << " , " 
					  << tRecoVtxPos[i].vtx.Z() << " , " 
					  << tRecoVtxPos[i].vtx.T() << ") " 
					  << "NLL = " << tRecoVtxPos[i].nll << std::endl;
		}
	}

	timer.Stop();
	fPerformances.vtx_search_ct = timer.RealTime();
	// std::cout << "SearchVertex took: " << timer.RealTime() << std::endl;

	double dStepSizeFinal = 10;
	int iToleranceFinal = 1;

	timer.Reset();
	timer.Start();
	
	ROOT::Math::XYZTVector limits(1.,1.,1.,1.); 
	limits *= 2.*LEAFConfig::fSearchVtxStep;
	fRecoVtxPosFinal = this->MinimizeVertex_Main(tRecoVtxPos, limits, dStepSizeFinal, LEAFConfig::fSearchVtxTolerance, iToleranceFinal, VERBOSE, true, false, LEAFConfig::fMinimizeLimitsNegative, LEAFConfig::fMinimizeLimitsPositive, LEAFConfig::fUseDirectionality, FilterThreshold);

	timer.Stop();
	fPerformances.vtx_minimize_ct = timer.RealTime();

	fOutput.vtx = fRecoVtxPosFinal[0].vtx;
	fOutput.vtx_nll = fRecoVtxPosFinal[0].nll;
}

std::vector<LEAFStructures::VtxCandidate> LEAF::SearchVertex_Main(int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality) {
// 1) Launch threads just like in SearchVertex_Main.
	std::vector<std::thread> lThreadList;
	std::vector<LEAFStructures::VtxCandidate> lOutputFinal;
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	int iCand_Step = (int)fPositionList.size()/fThread;
	if(iCand_Step < 1) iCand_Step = 1;

	LEAFConfig::mtx.lock();
	fVtxThreadOutput.clear();
	LEAFConfig::mtx.unlock();

	for(int iStart=0; iStart<(int)fPositionList.size(); iStart+=iCand_Step) {
		std::thread tThrd(
			LEAFThreads::SearchVertex_thread,
			iStart, iCand_Step,
			tolerance,
			likelihood, lowerLimit, upperLimit, directionality
		);
		lThreadList.push_back(std::move(tThrd));
	}

	// std::cout << "time after vtx search threads in SearchVertex_Main : " << timer.RealTime() << std::endl;

	for(auto &th: lThreadList) if(th.joinable()) th.join();

	// 2) collect partial results from fThreadOutput.
	LEAFConfig::mtx.lock();
	lOutputFinal = fVtxThreadOutput; 
	fVtxThreadOutput.clear();
	LEAFConfig::mtx.unlock();

	// 3) Sort by NLL and keep 'tolerance' number of candidates
	std::sort(lOutputFinal.begin(), lOutputFinal.end(), LEAF::SortingNLL<LEAFStructures::VtxCandidate>());
	while((int)lOutputFinal.size() > tolerance) lOutputFinal.pop_back();

	//4) Now build a vector of VtxCandidate.
	std::vector<LEAFStructures::VtxCandidate> finalList;
	finalList.reserve(lOutputFinal.size());

	for(auto &cand : lOutputFinal) {
		finalList.push_back(cand);
	}

	return finalList;
}

// Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
// The step between each grid points is given by stepSize (in centimeters).
// We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
// The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower).
// The number of output vertices can be chosen using "tolerance".
void LEAFThreads::SearchVertex_thread(int iStart, int iIte, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality) {
	LEAF::GetME()->SearchVertex_thread(iStart, iIte, tolerance, likelihood, lowerLimit, upperLimit, directionality);
}

void LEAF::SearchVertex_thread(int iStart, int iIte, int tolerance, bool likelihood, double lowerLimit, double upperLimit, int directionality) {
	std::vector<LEAFStructures::VtxCandidate> tVtxContainer;

	TStopwatch timer;
	timer.Reset();
	timer.Start();

	if (fPositionList.size() == 0)
		std::cout << " Position list not filled " << std::endl;

	int iEnd = (iIte + iStart);
	if (iEnd > (int)fPositionList.size()) iEnd = fPositionList.size();

	for (int iPos = iStart; iPos < iEnd; iPos++) {

		LEAFStructures::VtxCandidate hPos;
		hPos.vtx = fPositionList[iPos];
		hPos.nll = 0;

		if (likelihood) hPos.nll = this->Vertex_Time_NLL(hPos.vtx, lowerLimit, upperLimit, false, false, directionality);
		else hPos.nll = this->Vertex_Score(hPos.vtx, lowerLimit, upperLimit, false, false, directionality);

		//if (iPos < 10)
		//	std::cout << iPos << " -> " << hPos.vtx.X() << ", " <<  hPos.vtx.Y() << ", " <<  hPos.vtx.Z() << ", " <<  hPos.vtx.T() << " = " << hPos.nll << " " << likelihood << std::endl;

		tVtxContainer.push_back(hPos);
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), LEAF::SortingNLL<LEAFStructures::VtxCandidate>());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();

	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	LEAFConfig::mtx.lock();
	for (unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++) {
		fVtxThreadOutput.push_back(tVtxContainer[iPos]);
	}
	LEAFConfig::mtx.unlock();
}

std::vector<LEAFStructures::VtxCandidate> LEAF::MinimizeVertex_Main(std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality, float FilterThreshold) {
	std::vector<std::thread> lThreadList;
	std::vector<LEAFStructures::VtxCandidate> lOutputFinal;

	int iCand_Step = nCandidates / fThread;
	if (iCand_Step < 1) iCand_Step = 1;

	LEAFConfig::mtx.lock();
	fVtxThreadOutput.clear();
	LEAFConfig::mtx.unlock();

	for (int iStart = 0; iStart < nCandidates; iStart += iCand_Step) {
		std::thread tThrd(LEAFThreads::MinimizeVertex_thread, iStart, iCand_Step,
						  initialVertex, limits, stepSize, nCandidates,
						  tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality, FilterThreshold);

		lThreadList.push_back(std::move(tThrd));
	}

	for(auto &th: lThreadList) if(th.joinable()) th.join();

	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	LEAFConfig::mtx.lock();
	lOutputFinal = fVtxThreadOutput;
	fVtxThreadOutput.clear();
	LEAFConfig::mtx.unlock();

	std::sort(lOutputFinal.begin(), lOutputFinal.end(), LEAF::SortingNLL<LEAFStructures::VtxCandidate>());

	while ((int)lOutputFinal.size() > tolerance) lOutputFinal.pop_back();

	return lOutputFinal;
}

// here is the function where minimization gonna take place. The loop is inside this function
void LEAFThreads::MinimizeVertex_thread(int iStart, int iIte, std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates, int tolerance, int verbose, bool likelihood, bool average, double lowerLimit, double upperLimit, int directionality, float FilterThreshold) {
	LEAF::GetME()->MinimizeVertex_thread(iStart, iIte, initialVertex, limits, stepSize, nCandidates, tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality, FilterThreshold);
}

void LEAF::MinimizeVertex_thread(int iStart, int iIte, std::vector<LEAFStructures::VtxCandidate> initialVertex, ROOT::Math::XYZTVector limits, double stepSize, int nCandidates, int tolerance, int verbose, bool /*likelihood*/, bool /*average*/, double lowerLimit, double upperLimit, int directionality, float FilterThreshold) {
	int iEnd = (iIte + iStart);
	if (iEnd > nCandidates) iEnd = nCandidates;

	std::vector<LEAFStructures::VtxCandidate> tVtxContainer;

	LEAFConfig::mtx.lock();
	TFitter *minimizer = new TFitter(14); // 4=nb de params?
	LEAFConfig::mtx.unlock();
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
	minimizer->SetFCN(LEAFLikelihoods::MinuitLikelihood);		// here is function to minimize

	double Mig[2] = {1e6, 1e0}; // maxcalls and tolerance

	minimizer->SetParameter(0, "vertex0", 0, stepSize, -limits.X(), limits.X());
	minimizer->SetParameter(1, "vertex1", 0, stepSize, -limits.Y(), limits.Y());
	minimizer->SetParameter(2, "vertex2", 0, stepSize, -limits.Z(), limits.Z());
	minimizer->SetParameter(3, "vertex3", 0, stepSize / LEAFConfig::fLightSpeed, -limits.T() / LEAFConfig::fLightSpeed, limits.T() / LEAFConfig::fLightSpeed);

	minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
	minimizer->SetParameter(5, "nhits", 0, 0, -1, 1);
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

	for (int icand = iStart; icand < iEnd; icand++) {
		minimizer->SetParameter(0, "vertex0", initialVertex[icand].vtx.X(), stepSize, initialVertex[icand].vtx.X() - limits.X(), initialVertex[icand].vtx.X() + limits.X());
		minimizer->SetParameter(1, "vertex1", initialVertex[icand].vtx.Y(), stepSize, initialVertex[icand].vtx.Y() - limits.Y(), initialVertex[icand].vtx.Y() + limits.Y());
		minimizer->SetParameter(2, "vertex2", initialVertex[icand].vtx.Z(), stepSize, initialVertex[icand].vtx.Z() - limits.Z(), initialVertex[icand].vtx.Z() + limits.Z());
		minimizer->SetParameter(3, "vertex3", initialVertex[icand].vtx.Y(), stepSize / LEAFConfig::fLightSpeed, initialVertex[icand].vtx.T() - limits.T() / LEAFConfig::fLightSpeed, initialVertex[icand].vtx.T() + limits.T() / LEAFConfig::fLightSpeed);

		minimizer->SetParameter(4, "PMTConfiguration", 0, 0, 0, 0); // Useless now
		minimizer->SetParameter(5, "nhits", 0, 0, -1, 1);
		minimizer->SetParameter(6, "lowerLimit", lowerLimit, 5e-1, -7, -2);
		minimizer->SetParameter(7, "upperLimit", upperLimit, 5e-1, 2, 10);
		minimizer->SetParameter(8, "expoSigma", 100, 5, 0, 1e3);
		minimizer->SetParameter(9, "directionality", directionality, 1, 0, 2);

		if(fDirectionFilterSet) {
			minimizer->SetParameter(10, "dir_x", fDirectionFilter.X(), 5e-1, fDirectionFilter.X() - 1, fDirectionFilter.X() + 1);
			minimizer->SetParameter(11, "dir_y", fDirectionFilter.Y(), 5e-1, fDirectionFilter.Y() - 1, fDirectionFilter.Y() + 1);
			minimizer->SetParameter(12, "dir_z", fDirectionFilter.Z(), 5e-1, fDirectionFilter.Z() - 1, fDirectionFilter.Z() + 1);
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

		LEAFStructures::VtxCandidate hPos;
		hPos.vtx = ROOT::Math::XYZTVector(minimizer->GetParameter(0), minimizer->GetParameter(1), minimizer->GetParameter(2), minimizer->GetParameter(3));
		
		// Get minuit
		minuit = minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);

		hPos.nll = fmin;

		tVtxContainer.push_back(hPos);

		if (VERBOSE >= 2) {
			std::cout << "Initial vertex position = (" 
					  << initialVertex[icand].vtx.X() << "," 
					  << initialVertex[icand].vtx.Y() << "," 
					  << initialVertex[icand].vtx.Z() << "," 
					  << initialVertex[icand].vtx.T() << "), "
					  << "Candidate vertex position = (" 
					  << hPos.vtx.X() << "," 
					  << hPos.vtx.Y() << "," 
					  << hPos.vtx.Z() << "," 
					  << hPos.vtx.T() << "), "
					  << "NLL = " << hPos.nll << std::endl;
		}
	}

	std::sort(tVtxContainer.begin(), tVtxContainer.end(), LEAF::SortingNLL<LEAFStructures::VtxCandidate>());
	while ((int)tVtxContainer.size() > tolerance) tVtxContainer.pop_back();
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	LEAFConfig::mtx.lock();
	for (unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++) {
		fVtxThreadOutput.push_back(tVtxContainer[iPos]);
	}
	delete minimizer;
	LEAFConfig::mtx.unlock();
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
	
	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for(unsigned int i = 0; i < fHitsCollection[iPMTType]->size(); i++) {
			const HKHit* lHit = fHitsCollection[iPMTType]->at(i);
		
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
			
			totalPECorr += lHit->GetCharge();
			totalPE += lHit->GetCharge();
		}
	}
		
	double energyEstimate = totalPECorr * 0.111247 - 12.970173; //* Affine correction (affine because of dark rate and systematic errors), convert from charge to energy
	fOutput.total_charge = totalPE;
	fOutput.energy = energyEstimate;

	timer.Stop();
	fPerformances.energy_fit_ct = timer.RealTime();
}


/*****************************************************************************************************/
/* DIRECTION AND VERTEX FIT */
/*****************************************************************************************************/

std::vector<LEAF::JointFitCandidate> LEAF::SearchVertexAndDir() {
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	std::vector<LEAFStructures::VtxCandidate> tVtxCandidates = this->SearchVertex_Main(LEAFConfig::fSearchVtxTolerance,false,LEAFConfig::fSTimePDFLimitsQueueNegative,LEAFConfig::fSTimePDFLimitsQueuePositive,false);
	
	std::vector<JointFitCandidate> Candidates;
	for(int i=0; i<(int)tVtxCandidates.size(); i++) {
		std::vector<DirectionCandidate> tDirCandidates = this->SearchDirection(tVtxCandidates[i].vtx, LEAFConfig::fDirTolerance);
		for(int j=0; j<(int) tDirCandidates.size(); j++) {
			JointFitCandidate candidate;
			candidate.vtx = tVtxCandidates[i].vtx;
			candidate.theta = tDirCandidates[j].theta;
			candidate.phi = tDirCandidates[j].phi;
			candidate.nll = tDirCandidates[j].nll * tVtxCandidates[i].nll;
			Candidates.push_back(candidate);
		}
	}

	std::sort(Candidates.begin(), Candidates.end(), LEAF::SortingNLL<JointFitCandidate>{});

	while ((int)Candidates.size() > LEAFConfig::fSearchVtxTolerance) Candidates.pop_back();

	timer.Stop();
	fPerformances.vtx_search_ct = timer.RealTime();

	return Candidates;
}

LEAF::JointFitCandidate LEAF::MinimizeVertexAndDir(std::vector<JointFitCandidate> initialCandidates) {
	TStopwatch timer;
	timer.Reset();
	timer.Start();

	double stepSize = 10;
	ROOT::Math::XYZTVector limits(1.,1.,1.,1.); 
	limits *= 2.*LEAFConfig::fSearchVtxStep;

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
	minimizer->SetFCN(LEAFLikelihoods::MinuitJointNLL);		// here is function to minimize

	double Mig[2] = {1e6, 1e0}; // maxcalls and tolerance

	minimizer->SetParameter(0, "vertex0", 0, stepSize, -limits.X(), limits.X());
	minimizer->SetParameter(1, "vertex1", 0, stepSize, -limits.Y(), limits.Y());
	minimizer->SetParameter(2, "vertex2", 0, stepSize, -limits.Z(), limits.Z());
	minimizer->SetParameter(3, "vertex3", 0, stepSize / LEAFConfig::fLightSpeed, -limits.T() / LEAFConfig::fLightSpeed, limits.T() / LEAFConfig::fLightSpeed);
	minimizer->SetParameter(4, "theta", 0, 2*M_PI/180, 0.0, TMath::Pi());
	minimizer->SetParameter(5, "phi", 0, 2*M_PI/180, -TMath::Pi(), TMath::Pi());

	minimizer->SetParameter(6, "PMTConfiguration", 0, 0, 0, 0); // Useless now
	minimizer->SetParameter(7, "nhits", 0, 0, -1, 1);
	minimizer->SetParameter(8, "lowerLimit", LEAFConfig::fSTimePDFLimitsQueueNegative, 5e-1, -7, -2);
	minimizer->SetParameter(9, "upperLimit", LEAFConfig::fSTimePDFLimitsQueuePositive, 5e-1, 2, 10);
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
		minimizer->SetParameter(0, "vertex0", initialCandidates[icand].vtx.X(), stepSize, initialCandidates[icand].vtx.X() - limits.X(), initialCandidates[icand].vtx.X() + limits.X());
		minimizer->SetParameter(1, "vertex1", initialCandidates[icand].vtx.Y(), stepSize, initialCandidates[icand].vtx.Y() - limits.Y(), initialCandidates[icand].vtx.Y() + limits.Y());
		minimizer->SetParameter(2, "vertex2", initialCandidates[icand].vtx.Z(), stepSize, initialCandidates[icand].vtx.Z() - limits.Z(), initialCandidates[icand].vtx.Z() + limits.Z());
		minimizer->SetParameter(3, "vertex3", initialCandidates[icand].vtx.T(), stepSize / LEAFConfig::fLightSpeed, initialCandidates[icand].vtx.T() - limits.T() / LEAFConfig::fLightSpeed, initialCandidates[icand].vtx.T() + limits.T() / LEAFConfig::fLightSpeed);
		minimizer->SetParameter(4, "theta", initialCandidates[icand].theta, 2*M_PI/180, 0.0, TMath::Pi());
		minimizer->SetParameter(5, "phi", initialCandidates[icand].phi, 2*M_PI/180, -TMath::Pi(), TMath::Pi());

		minimizer->ExecuteCommand("MIGRAD", Mig, 2);

		JointFitCandidate tCandReco;

		tCandReco.vtx.SetXYZT(minimizer->GetParameter(0), minimizer->GetParameter(1), minimizer->GetParameter(2), minimizer->GetParameter(3));
		tCandReco.theta = minimizer->GetParameter(4);
		tCandReco.phi = minimizer->GetParameter(5);

		// Get minuit
		double fmin, fedm, errdef;
		int npar, nparx, istat;
		minuit->mnstat(fmin, fedm, errdef, npar, nparx, istat);
		tCandReco.nll = fmin;

		tContainer.push_back(tCandReco);
	}

	std::sort(tContainer.begin(), tContainer.end(), LEAF::SortingNLL<JointFitCandidate>{});

	delete minimizer;

	timer.Stop();
	fPerformances.vtx_minimize_ct = timer.RealTime();

	return tContainer[0];
}

