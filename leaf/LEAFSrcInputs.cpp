/*****************************************************************************************************/
/**	LEAFSrcInputs.cc																				**/
/**	Date: March 27th 2026																			**/
/**	Desc: Implementation of the I/O functions													**/
/*****************************************************************************************************/

#ifndef LEAF_HPP
#include "LEAF.hpp"
#endif

#ifdef HK_USE_ROOT7
void LEAF::Initialize(std::shared_ptr<HKGeometry const> lGeometry, std::shared_ptr<HKDarkNoise const> lDarkNoise, std::shared_ptr<HKGeometryPMTCollection const> lGeoPMT_ID, std::shared_ptr<HKGeometryPMTCollection const> lGeoPMT_mPMT)
#else
void LEAF::Initialize(const HKGeometry *lGeometry, const HKDarkNoise *lDarkNoise, const HKGeometryPMTCollection *lGeoPMT_ID, const HKGeometryPMTCollection *lGeoPMT_mPMT)
#endif
{
    fThread = N_THREAD;
    fRand = new TRandom3();
	
	fGeometry = lGeometry;
	fDarkNoise = lDarkNoise;
	fGeoPMTs[int(PMTType::kID)] = lGeoPMT_ID;
	fGeoPMTs[int(PMTType::kmPMT)] = lGeoPMT_mPMT;

	LEAFConfig::Initialize(fGeometry->GetInnerTyvekRadius(),fGeometry->GetInnerTyvekTotalHeight());
	this->InitInputs();
	this->LoadSplines();
}

void LEAF::InitInputs() {
    this->MakePositionList();
}

void LEAF::LoadSplines() {
	std::cout << "Loading splines..." << std::endl;
	TFile *fSplines;
	if (LEAFConfig::fHighEnergy) {
		fSplines = new TFile("${LEAFDIR}/inputs/PDF_electron_upto100MeV_uniform_isotropic_10k_withDN.root", "read");
	}
	else {
		fSplines = new TFile("${LEAFDIR}/inputs/PDF_electron_upto100MeV_uniform_isotropic_10k_withDN.root","read");
	}

	// Prevent TGraph2D to be append to TFile (this is needed as we are doing multiple copy of TGraph2D)
	TH1::AddDirectory(false);

	//Check if the files are opened
	if (!fSplines || !fSplines->IsOpen()) {
		std::cerr << "Error: Splines file is not opened." << std::endl;
		exit(1);
	}

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		int iPMTType = static_cast<int>(pmtType);

		std::cout << "Process spline" << std::endl;

		//Load TOF residuals splines
		fSplineTimePDFConv[iPMTType]     = (TSpline3*) fSplines->Get(Form("splineExpoConv_%d", iPMTType));
		fSplineTimePDFQueue[iPMTType]    = (TSpline3*) fSplines->Get(Form("splineExpoQueue_%d", iPMTType));
		fSplineTimePDFDarkRate[iPMTType] = (TSpline3*) fSplines->Get(Form("splineDR_%d", iPMTType));
		fDirectionPDF[iPMTType]          = (TSpline3*) fSplines->Get(Form("splineChargeProfile_%d", iPMTType));

		// std::cout << "ok ?" << std::endl;
		LEAFConfig::fSTimePDFLimitsQueueNegative_fullTimeWindow = fSplineTimePDFQueue[iPMTType]->GetXmin();
		LEAFConfig::fSTimePDFLimitsQueuePositive_fullTimeWindow = fSplineTimePDFQueue[iPMTType]->GetXmax();
		//std::cout<<"Min="<<fSTimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<fSTimePDFLimitsQueuePositive_fullTimeWindow<<std::endl;
		if(LEAFConfig::fSTimePDFLimitsQueueNegative < LEAFConfig::fSTimePDFLimitsQueueNegative_fullTimeWindow) LEAFConfig::fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[iPMTType]->GetXmin();//To avoid to use PDF where it is not defined
		if(LEAFConfig::fSTimePDFLimitsQueuePositive > LEAFConfig::fSTimePDFLimitsQueuePositive_fullTimeWindow) LEAFConfig::fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[iPMTType]->GetXmax();//same here.

		// std::cout << "okay ?" << std::endl;
		//Load 3D directionality histograms.
		std::cout << "directionality" << std::endl;
	}
	std::cout << "Splines loaded." << std::endl;
}

void LEAF::MakePositionList() {
	// Make position list for a given step size
	double dStep = LEAFConfig::fSearchVtxStep;
	fPositionList.clear();

	for (double time = -50; time < 50; time += (dStep / LEAFConfig::fLightSpeed)) {
		for (double radius = 0; radius <= LEAFConfig::fTankRadius; radius += dStep) {
			double perimeter = 2 * TMath::Pi() * radius;
			int numberOfStepsOnCircle = floor(perimeter / dStep);
			if (numberOfStepsOnCircle == 0) numberOfStepsOnCircle = 1;
			double angleStepOnCircle = 2 * TMath::Pi() / numberOfStepsOnCircle;

			for (double angle = 0; angle <= 2 * TMath::Pi(); angle += angleStepOnCircle) {
				for (double height = -LEAFConfig::fTankHalfHeight; height <= LEAFConfig::fTankHalfHeight; height += dStep) {
					ROOT::Math::XYZTVector vertex(
						radius * TMath::Cos(angle),
						radius * TMath::Sin(angle),
						height,
						time
					); // In centimeters
					fPositionList.push_back(vertex);
				}
			}
		}
	}
}

#ifdef HK_USE_ROOT7
void LEAF::LoadHitsCollection(std::shared_ptr<HKHitsCollection const> lHitCol_ID, std::shared_ptr<HKHitsCollection const> lHitCol_mPMT) 
#else
void LEAF::LoadHitsCollection(const HKHitsCollection *lHitCol_ID, const HKHitsCollection * lHitCol_mPMT) 
#endif
{
	fHitsCollection[int(PMTType::kID)] = lHitCol_ID;
	fHitsCollection[int(PMTType::kmPMT)] = lHitCol_mPMT;

	if (fPositionList.size() == 0) MakePositionList();
	if (fHitsCollection[int(PMTType::kmPMT)] == nullptr) {
		LEAFConfig::fUseDirectionality = false;
	}
}

void LEAF::SetTrueVertexInfo(ROOT::Math::XYZTVector vtx) {
	fTrueVtxPos = vtx;
}

void LEAF::SetTrueDirInfo(ROOT::Math::XYZVector  trueDir) {
	fTrueDir = trueDir;
}

void LEAF::SetNThread(int iThread) {
	fThread = iThread; 
}

void LEAF::GenerateEventInfo(double lowerLimit, double upperLimit) {
	
	// Check if process was already completed for these limits
	if (fLastLowerLimit == lowerLimit && fLastUpperLimit == upperLimit) return;

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		int iPMTType = static_cast<int>(pmtType);

		fLastLowerLimit = lowerLimit;
		fLastUpperLimit = upperLimit;

		// Initialize
		fEventInfo[iPMTType].nHits = fHitsCollection[iPMTType]->size();
		fEventInfo[iPMTType].SignalOverNoise = 0.;
		fEventInfo[iPMTType].NoiseIntegral = 0.;
		fEventInfo[iPMTType].SignalIntegral = 0.;

		double signalDR, signalPE;
		signalDR = LEAFConfig::fTimeWindowSizeFull * fDarkRate_ns[iPMTType];		// Over the whole time window
		signalPE = std::max(fEventInfo[iPMTType].nHits - signalDR, 0.); // Over the whole time window.

		double signalPETime = signalPE; // I assume that all signal is here.
		double signalDRTime = (upperLimit - lowerLimit) * signalDR / LEAFConfig::fTimeWindowSizeFull;

		fEventInfo[iPMTType].SignalOverNoise = signalPETime / signalDRTime; // Over the
		fEventInfo[iPMTType].NoiseIntegral = SplineUtilities::SplineIntegral(fSplineTimePDFDarkRate[iPMTType], lowerLimit, upperLimit);
		fEventInfo[iPMTType].SignalIntegral = SplineUtilities::SplineIntegralAndSubstract(fSplineTimePDFQueue[iPMTType], fSplineTimePDFDarkRate[iPMTType], lowerLimit, upperLimit);

		if (VERBOSE >= 3) {
			std::cout << "nhits=" << fEventInfo[iPMTType].nHits 
					  << ", DR average=" << signalDR
					  << ", signal over noise=" << fEventInfo[iPMTType].SignalOverNoise
					  << ", signal integral=" << fEventInfo[iPMTType].SignalIntegral
					  << ", DR integral=" << fEventInfo[iPMTType].NoiseIntegral
					  << ", in integral=" << fEventInfo[iPMTType].SignalIntegral / fEventInfo[iPMTType].NoiseIntegral 
					  << std::endl;
			
		} 
	}
}

double LEAF::CalculateSNR(ROOT::Math::XYZTVector vertex, double lowerLimit, double upperLimit) {

	// We'll do a simple approach: count how many hits are in [lower,upper]
	// for each PMT type, then subtract the dark rate.

	int nHitsTotal = 0;
	std::array<int, static_cast<size_t>(PMTType::kNumPMTTypes)> inTimes;

	ROOT::Math::XYZVector vertexPos(vertex.X(),vertex.Y(),vertex.Z());
	
	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		nHitsTotal += fHitsCollection[iPMTType]->size();

		int nHitsInTime = 0;
		for(unsigned int i = 0; i < fHitsCollection[iPMTType]->size(); i++) {
			std::vector<double> toHit;
			
			const HKHit* lHit = fHitsCollection[iPMTType]->at(i);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);
	
			double hitTime = lHit->GetTime();
			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, hitTime, vertex.T());

			//std::cout << "residual: " << residual << ", " << "lowerlimit and upperlimit " << lowerLimit << ", " << upperLimit << std::endl; 
			if(residual >= lowerLimit && residual <= upperLimit) {
				nHitsInTime += 1;
			}
		}
		inTimes[iPMTType] = nHitsInTime;
	}
	
	double dt = (upperLimit - lowerLimit);
	
	// BnL:
	double darkRateBnL = fDarkNoise->GetAverageTotalDarkNoisePerNS(PMTType::kID);
	// average dark over dt:
	double darkInWindowBnL = dt * darkRateBnL;
	double signalBnL = std::max(double(inTimes[int(PMTType::kID)]) - darkInWindowBnL, 0.0);
	//std::cout << "inTimeBnL - darkInWindowBnL : " << inTimeBnL << "-" << darkInWindowBnL << std::endl;
	
	double snr = (darkInWindowBnL > 1e-9) ? (signalBnL / darkInWindowBnL) : 999.0;
	return snr;
}