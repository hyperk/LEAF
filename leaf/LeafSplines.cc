#include "LeafSplines.hh"

TSpline3 *	fSplineTimePDFQueue[NPMT_CONFIGURATION];
TSpline3 *	fSplineTimePDFDarkRate[NPMT_CONFIGURATION];
TSpline3* fDirectionPDF;
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
	TFile *fSplines, *fSplines2, *fSplines1;
	if (fHighEnergy)
	{
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_HE.root", "read"); // To generate with my code ProduceWSPlots.c
		fSplines2 = fSplines;
	}
	else
	{
		std::cout << "spline 1" << std::endl;
		fSplines1 = new TFile("${LEAFDIR}/inputs/timePDFNoDR_50000_e10MeV_Hit_720.root","read");//To generate with my code ProduceWSPlots.c		
		// fSplines = new TFile("${LEAFDIR}/inputs/timePDF_DRnew_Large.root", "read");			  // To generate with my code ProduceWSPlots.c
		std::cout << "spline 2" << std::endl;
		// fSplines = new TFile("${LEAFDIR}/inputs/timePDFDR_5000_e10MeV_Hit_fiducial.root","read"); //*TAHA
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_3M_10T_500000.root","read"); //*NICOLAS
		// fSplines = new TFile("${LEAFDIR}/inputs/timePDF_3M_10T.root","read");
		std::cout << "spline 3" << std::endl;
		fSplines2 = new TFile("${LEAFDIR}/inputs/timePDF_Directionality_DRnew.root", "read"); // To generate with my code ProduceWSPlots.c
		// fSplines2 = new TFile("${LEAFDIR}/inputs/timePDF_Directionality_DRnew.root","read");//To generate with my code ProduceWSPlots.c
	}

	std::cout << "spline file read" << std::endl;

	// Prevent TGraph2D to be append to TFile (this is needed as we are doing multiple copy of TGraph2D)
	// For a strange reason sometimes the TGraph2D is see as a TH1 and stored in an open TFile
	TH1::AddDirectory(false);

	//Check if the files are opened
	if (!fSplines || !fSplines->IsOpen())
	{
		std::cerr << "Error: Splines file is not opened." << std::endl;
		exit(1);
	}
	if (!fSplines2 || !fSplines2->IsOpen())
	{
		std::cerr << "Error: Splines file is not opened." << std::endl;
		exit(1);
	}
	if (!fSplines1 || !fSplines1->IsOpen())
	{
		std::cerr << "Error: Splines file is not opened." << std::endl;
		exit(1);
	}

	// std::cout << "splines files loaded" << std::endl;

	int configs = NPMT_CONFIGURATION;

	for(int pmtType=0;pmtType<configs;pmtType++)
	{
		std::cout << "Get spline" << std::endl;
		//Load 1D t-tof splines
		fSplineTimePDFQueue[pmtType]    = (TSpline3*) fSplines->Get(Form("splineExpoQueue%d_%d",0,pmtType));
		fSplineTimePDFDarkRate[pmtType] = (TSpline3*) fSplines->Get(Form("splineDR%d_%d",0,pmtType));
		fDirectionPDF = (TSpline3*) fSplines1->Get(Form("ChargeProfileCos_PDF_spline_pmtType0"));
		if (!fDirectionPDF) 
		{
			std::cerr << "Error: fDirectionPDF is not loaded properly." << std::endl;
			exit(1);
		}

		std::cout << "process spline" << std::endl;
		// std::cout << "ok ?" << std::endl;
		fSTimePDFLimitsQueueNegative_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmin();
		fSTimePDFLimitsQueuePositive_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmax();
		//std::cout<<"Min="<<fSTimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<fSTimePDFLimitsQueuePositive_fullTimeWindow<<std::endl;
		if(fSTimePDFLimitsQueueNegative < fSTimePDFLimitsQueueNegative_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmin();//To avoid to use PDF where it is not defined
		if(fSTimePDFLimitsQueuePositive > fSTimePDFLimitsQueuePositive_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmax();//same here.

		// std::cout << "okay ?" << std::endl;
		//Load 3D directionality histograms.
		std::cout << "directionality" << std::endl;
		// Ensure HKAA::kmPMT_Groups is defined and initialized
		if (!fSplines2 || !fSplines2->IsOpen()) {
			std::cerr << "Error: fSplines2 is not opened or is null." << std::endl;
			exit(1);
		}

		int grpNB = 0;
		if (HKAA::kmPMT_Groups > 0) {
			grpNB = HKAA::kmPMT_Groups;
		} else {
			std::cerr << "Error: HKAA::kmPMT_Groups is not properly defined or initialized." << std::endl;
			exit(1);
		}

		std::cout << "got grnNB" << std::endl;

		for (int pmtGroup = 0; pmtGroup < grpNB; pmtGroup++) {
			gPMTDirectionality_2D[pmtType][pmtGroup] = (TGraph2D*) fSplines2->Get(Form("gPMTDirectionality_2D_%d_%d_%d", 0, pmtType, pmtGroup));
			if (!gPMTDirectionality_2D[pmtType][pmtGroup]) {
				std::cerr << "Error: Failed to load gPMTDirectionality_2D for pmtType " << pmtType << " and pmtGroup " << pmtGroup << std::endl;
				exit(1);
			}
		}
		
		std::cout << "dist response" << std::endl;

		fDistResponsePMT[pmtType] = (TF1*) fSplines2->Get(Form("fDistResponsePMT_pmtType%d", pmtType));
		if (!fDistResponsePMT[pmtType]) 
		{
			std::cerr << "Error: Failed to load fDistResponsePMT for pmtType " << pmtType << std::endl;
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