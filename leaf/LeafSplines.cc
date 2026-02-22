#include "LeafSplines.hh"

TSpline3 *	fSplineTimePDFConv[NPMT_CONFIGURATION];
TSpline3 *	fSplineTimePDFQueue[NPMT_CONFIGURATION];
TSpline3 *	fSplineTimePDFDarkRate[NPMT_CONFIGURATION];
TSpline3 *  fDirectionPDF[NPMT_CONFIGURATION];
TGraph2D * 	gPMTDirectionality_2D[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];
TF1 * 		fDistResponsePMT[NPMT_CONFIGURATION];

double fDarkRate_dir_proba[NPMT_CONFIGURATION][HKAA::kmPMT_Groups];
double fLastLowerLimit;
double fLastUpperLimit;

double SplineIntegral(TSpline3 *s, double start, double end, double stepSize)
{
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * s->Eval(i + stepSize / 2);
	}
	return integral;
}
double SplineIntegralAndSubstract(TSpline3 *s0, TSpline3 *s1, double start, double end, double stepSize)
{
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * (s0->Eval(i + stepSize / 2) - s1->Eval(i + stepSize / 2));
	}
	return integral;
}

double SplineIntegralExpo(TSpline3 *s, double start, double end, double sigma, double stepSize)
{
	double integral = 0;
	for (double i = start; i < end; i += stepSize)
	{
		integral += stepSize * s->Eval(i + stepSize / 2) * TMath::Gaus(i + stepSize / 2, 0, sigma);
	}
	return integral;
}

void LoadSplines()
{
	std::cout << "Loading splines..." << std::endl;
	TFile *fSplines;
	if (fHighEnergy)
	{
		fSplines = new TFile("${LEAFDIR}/inputs/PDF_electron_upto100MeV_uniform_isotropic_10k_withDN.root", "read");
	}
	else
	{
		fSplines = new TFile("${LEAFDIR}/inputs/PDF_electron_upto100MeV_uniform_isotropic_10k_withDN.root","read");
	}

	// Prevent TGraph2D to be append to TFile (this is needed as we are doing multiple copy of TGraph2D)
	TH1::AddDirectory(false);

	//Check if the files are opened
	if (!fSplines || !fSplines->IsOpen())
	{
		std::cerr << "Error: Splines file is not opened." << std::endl;
		exit(1);
	}

	int configs = NPMT_CONFIGURATION;
	for(int pmtType=0 ; pmtType<configs ; pmtType++)
	{
		std::cout << "Process spline" << std::endl;

		//Load TOF residuals splines
		fSplineTimePDFConv[pmtType]     = (TSpline3*) fSplines->Get(Form("splineExpoConv_%d", pmtType));
		fSplineTimePDFQueue[pmtType]    = (TSpline3*) fSplines->Get(Form("splineExpoQueue_%d", pmtType));
		fSplineTimePDFDarkRate[pmtType] = (TSpline3*) fSplines->Get(Form("splineDR_%d", pmtType));
		fDirectionPDF[pmtType]          = (TSpline3*) fSplines->Get(Form("splineChargeProfile_%d", pmtType));

		// std::cout << "ok ?" << std::endl;
		fSTimePDFLimitsQueueNegative_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmin();
		fSTimePDFLimitsQueuePositive_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmax();
		//std::cout<<"Min="<<fSTimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<fSTimePDFLimitsQueuePositive_fullTimeWindow<<std::endl;
		if(fSTimePDFLimitsQueueNegative < fSTimePDFLimitsQueueNegative_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmin();//To avoid to use PDF where it is not defined
		if(fSTimePDFLimitsQueuePositive > fSTimePDFLimitsQueuePositive_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmax();//same here.

		// std::cout << "okay ?" << std::endl;
		//Load 3D directionality histograms.
		std::cout << "directionality" << std::endl;

		int grpNB = 0;
		if (HKAA::kmPMT_Groups > 0) 
		{
			grpNB = HKAA::kmPMT_Groups;
		}
		else
		{
			std::cerr << "Error: HKAA::kmPMT_Groups is not properly defined or initialized." << std::endl;
			exit(1);
		}
	}
	std::cout << "Splines loaded." << std::endl;
}

EventInfo MakeEventInfo(double lowerLimit, double upperLimit, int pmtType)
{
    struct EventInfo fEventInfo;

	if (fLastLowerLimit == lowerLimit && fLastUpperLimit == upperLimit) return fEventInfo;
	fLastLowerLimit = lowerLimit;
	fLastUpperLimit = upperLimit;

    fEventInfo.hits = 0;
    fEventInfo.SignaloverNoise = 0.;
    fEventInfo.NoiseIntegral = 0.;
    fEventInfo.SignalIntegral = 0.;

	int iHitTotal = fHitCollection->Size();

	for (int iHit = 0; iHit < iHitTotal; iHit++)
	{
		// Hit lHit = fHitInfo[iHit];
		Hit lHit = fHitCollection->At(iHit);

		int iPMT = lHit.PMT;
		int iType = Astro_GetPMTType(iPMT);
        if(pmtType == iType) fEventInfo.hits += 1;
	}

    double signalDR, signalPE;
    signalDR = fTimeWindowSizeFull * fDarkRate_ns[pmtType];		// Over the whole time window
    signalPE = std::max(fEventInfo.hits - signalDR, 0.); // Over the whole time window.

    double signalPETime = signalPE; // I assume that all signal is here.
    double signalDRTime = (upperLimit - lowerLimit) * signalDR / fTimeWindowSizeFull;

    fEventInfo.SignaloverNoise = signalPETime / signalDRTime; // Over the
    fEventInfo.NoiseIntegral = SplineIntegral(fSplineTimePDFDarkRate[pmtType], lowerLimit, upperLimit);
    fEventInfo.SignalIntegral = SplineIntegralAndSubstract(fSplineTimePDFQueue[pmtType], fSplineTimePDFDarkRate[pmtType], lowerLimit, upperLimit);

    if (VERBOSE >= 3) std::cout 
                << "nhits=" << fEventInfo.hits 
                << ", DR average=" << signalDR
                << ", signal over noise=" << fEventInfo.SignaloverNoise
                << ", signal integral=" << fEventInfo.SignalIntegral
                << ", DR integral=" << fEventInfo.NoiseIntegral
                << ",in integral=" << fEventInfo.SignalIntegral / fEventInfo.NoiseIntegral << std::endl;
        
    return fEventInfo;
}