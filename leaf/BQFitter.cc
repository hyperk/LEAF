/*********************************************************************************/
/**	BQFitter.cc								**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		**/
/**	Original author: Benjamin Quilain					**/
/**	Date: December 18th 2019						**/
/**	Desc: Low-E Fitter for Hyper-K						**/
/*********************************************************************************/

#include "BQFitter.hh"
#include "TStopwatch.h"

double BQFitter::fSTimePDFLimitsQueueNegative = 0;
double BQFitter::fSTimePDFLimitsQueuePositive = 0;
BQFitter* BQFitter::myFitter=NULL;

/************************************************************************************************************************/

void MinuitLikelihood(int& /*nDim*/, double * /*gout*/, double & NLL, double par[], int /*flg*/){
	//TStopwatch timer;
	//timer.Reset();
	//timer.Start();
	
	std::vector<double> vertexPosition(4,0.);//In centimeters
	for(int i=0;i<4;i++) vertexPosition[i]=par[i];
	//int pmtType=par[4];
	int nhits=par[5];
	double lowerLimit=par[6];
	double upperLimit=par[7];
	//double expoSigma=par[8];
	double directionality=par[9];
	//double scalingFactor=par[8];
	//NLL=findNLLSignalXDR(pmtType,vertexPosition,nhits,true,1,lowerLimit,upperLimit);
	//NLL=findNLLExpo(pmtType,vertexPosition,nhits,true,1,lowerLimit,upperLimit);
	//NLL=BQFitter::GetME()->FindNLL(vertexPosition,nhits,true,1,lowerLimit,upperLimit,true,false,directionality);
	NLL=BQFitter::GetME()->FindNLL_Likelihood(vertexPosition,nhits,lowerLimit,upperLimit,true,false,directionality);
	//cout<<"Directionality"<<endl;
	//NLL=findRateAndShape(pmtType,vertexPosition,nhits,1,lowerLimit,upperLimit,scalingFactor,true);
	//NLL=findChi2Binned(pmtType,vertexPosition,nhits,1,lowerLimit,upperLimit);
	//NLL=findChi2BinnedExpo(pmtType,vertexPosition,nhits,1,expoSigma);
	//if(verbose){
	//	cout<<setprecision(10)<<"NLL="<<NLL<<", Vertex: "<<vertexPosition[0]<<"ns, ("<<vertexPosition[1]<<", "<<vertexPosition[2]<<", "<<vertexPosition[3]<<")"<<endl;
	//	cout<<"Expo sigma = "<<expoSigma<<"ns"<<endl;
	//}
	//cout<<"Timing window = "<<lowerLimit<<"ns - "<<upperLimit<<"ns"<<endl<<endl;
	//double nll=BQFitter::GetME()->FindNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit);
	
	
	//timer.Stop();7
	//std::cout << " Likelihood: " << timer.RealTime() << std::endl;

}


void MinimizeVertex_CallThread(	
		int iStart, int iIte,	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
	
	BQFitter::GetME()->MinimizeVertex_thread(iStart, iIte,	
		initialVertex, limits, stepSize, nhits,
		nCandidates, tolerance, verbose, likelihood, average, lowerLimit,  upperLimit, directionality);
	
}	

void SearchVertex_CallThread(
			int iStart, int iIte,
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
			
	BQFitter::GetME()->SearchVertex_thread(iStart, iIte,	
		nhits, tolerance, likelihood, lowerLimit, upperLimit, directionality);
			
}
/************************************************************************************************************************/

BQFitter::BQFitter() {

	myFitter = this;
	this->Init();

}

BQFitter::~BQFitter() {


	//fPMT_TubeInMPMT.clear();
	//fPMT_Group.clear();
	//fPMT_Info.clear();
	fHitInfo.clear();
	
	fTrueVtxPosDouble.clear();
	delete fRand;
	
}

BQFitter* BQFitter::GetME() {

	if ( myFitter ) return myFitter;

	myFitter = new BQFitter();
	return myFitter;
}

/************************************************************************************************************************/

void BQFitter::Init() {

	// Parameters
	fStepByStep 					= false;
	fUseDirectionality 				= true;
	fHighEnergy					= false;
	
	// If true, we apply same cuts for timing as B&L PMT when searching vertex. 
	// Otherwise, we do not apply them and use all the hits. 
	// The latter is particularly useful when directionality is added, as it can kill DR hits.
	fLimit_mPMT 					= true; 
	
	fSTimePDFLimitsQueueNegative 			= -3.; //-5;//-3;//-1.5;//-3;//-10;
	fSTimePDFLimitsQueuePositive 			= 4.; //10;//4;//3;//4;//500;
	fSTimePDFLimitsQueueNegative_fullTimeWindow 	= 0.; //-1.5;//-3;//-10;
	fSTimePDFLimitsQueuePositive_fullTimeWindow 	= 0.; //3;//4;//500;
	fTimeWindowSizeFull 				= 1500;
	fAveraging 					= 20;//Number of points we used to average the vertex on.
	
	fMinimizeLimitsNegative				= -700;
	fMinimizeLimitsPositive				= 1000;
	
	fSearchVtxStep					= 300;
	
	fHitTimeLimitsNegative				= -5;//-10;//-5;//-5
	fHitTimeLimitsPositive				= 7;//15;//7;//7
	
	fIntegrationTimeWindow 				= 50;//in ns
	
	// Compute light speed:
	float fCVacuum 					= 3e8*1e2 / 1e9;//speed of light, in centimeter per ns.
	float fNIndex 					= 1.373;//1.385;//1.373;//refraction index of water
	fLightSpeed 					= fCVacuum/fNIndex;
	
	fTrueVtxPos.clear();
	fTrueVtxPosDouble.clear();
	fTrueVtxPos.resize(5,0.);
	fTrueVtxPosDouble.push_back(fTrueVtxPos);
	
	fPDFNorm_fullTimeWindow 			= 0.;
	
	//fPositionSize					= 0;
			
	// Get random generator
	fRand  						= new TRandom3();
	//fRand2 					= new TRandom3();
	
	this->LoadSplines();

}



void BQFitter::LoadSplines() {
	TFile *fSplines, *fSplines2;
	if(fHighEnergy){
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_HE.root","read");//To generate with my code ProduceWSPlots.c
		fSplines2 = fSplines;
	}
	else{
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_DRnew_Large.root","read");//To generate with my code ProduceWSPlots.c
		fSplines2 = new TFile("${LEAFDIR}/inputs/timePDF_Directionality.root","read");//To generate with my code ProduceWSPlots.c
	}
	for(int pmtType=0;pmtType<NPMT_CONFIGURATION;pmtType++){
		//Load 1D t-tof splines
		fSplineTimePDFQueue[pmtType]    = (TSpline3*) fSplines->Get(Form("splineExpoQueue%d_%d",0,pmtType));
		fSplineTimePDFDarkRate[pmtType] = (TSpline3*) fSplines->Get(Form("splineDR%d_%d",0,pmtType));

		fSTimePDFLimitsQueueNegative_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmin();
		fSTimePDFLimitsQueuePositive_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmax();
		std::cout<<"Min="<<fSTimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<fSTimePDFLimitsQueuePositive_fullTimeWindow<<std::endl;
		if(fSTimePDFLimitsQueueNegative < fSTimePDFLimitsQueueNegative_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmin();//To avoid to use PDF where it is not defined
		if(fSTimePDFLimitsQueuePositive > fSTimePDFLimitsQueuePositive_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmax();//same here.

		//Load 3D directionality histograms.
		for(int pmtGroup=0;pmtGroup<NGROUP_PMT;pmtGroup++){
			hPMTDirectionality_1D[pmtType][pmtGroup] = (TH1D*) fSplines2->Get(Form("hPMTDirectionality_1D_%d_%d_%d",0,pmtType,pmtGroup));
			fDarkRate_dir_proba[pmtType][pmtGroup] = hPMTDirectionality_1D[pmtType][pmtGroup]->Integral() / hPMTDirectionality_1D[pmtType][pmtGroup]->GetNbinsX();
			std::cout << hPMTDirectionality_1D[pmtType][pmtGroup]->GetMean() << ", and DR = "<<fDarkRate_dir_proba[pmtType][pmtGroup] << std::endl;
		}
		/*
		bohPMTDirectionality_1Dol testSplines=true;
		if(testSplines){
			std::cout<<"testing spline for pmt type = " << pmtType << std::endl;
			std::cout<<"Min value of the spline range = " << stimePDF[pmtType]->GetXmin() << std::endl;
			for(double t=-stimePDFLimits;t<stimePDFLimits;t+=1){
				std::cout << "Time " << t << "ns, value = "<< stimePDF[pmtType]->Eval(t) << std::endl;
				//std::cout << "Time " << t << "ns, value = "<< ftimePDF[pmtType]->Eval(t) << std::endl;
			}
		}
		ftimePDF[pmtType] = new TF1(Form("ftimePDF%d",pmtType),valueTimePDF,-ftimePDFLimits,ftimePDFLimits,3);
		ftimePDF[pmtType]->SetParameter(0,-20);
		ftimePDF[pmtType]->SetParameter(1,20);
		ftimePDF[pmtType]->SetParameter(2,pmtType);
		*/
	}
}

void BQFitter::SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal, double dDarkRate_mPMT ) {

	// Keep RootGeometry event:
	fWCGeo = wGeo;
	
	// Read important variables:
	fTankRadius = fWCGeo->GetWCCylRadius(); //Conversion from cm to m removed
	fTankHeight = fWCGeo->GetWCCylLength(); //Conversion from cm to m removed
	fTankHalfHeight = fTankHeight / 2.;
	
	fDarkRate_Normal = dDarkRate_Normal;
	fDarkRate_mPMT   = dDarkRate_mPMT;
	
	if ( dDarkRate_Normal == 0 ) {
		fDarkRate_ns[0] = 8400.  * 2e4 * /*fWCGeo->GetWCNumPMT(false) */ 1e-9;
		fDarkRate_ns[1] = 100.   * 5e3 * /*fWCGeo->GetWCNumPMT(true ) */ 1e-9 * 19;
	}
	else {
		fDarkRate_ns[0] = fDarkRate_Normal * fWCGeo->GetWCNumPMT(false) * 1e-9;
		fDarkRate_ns[1] = fDarkRate_mPMT   * fWCGeo->GetWCNumPMT(true ) * 1e-9 * 19;
	}
	
	this->LoadPMTInfo();
	this->LoadDirConstant();
	
	// Make position lists
      	this->MakePositionList();
}

void BQFitter::SetTrueVertexInfo(std::vector<double> vtx, double time) {

	fTrueVtxPosDouble.clear();
	fTrueVtxPos.clear();
	
	fTrueVtxPos.resize(5,0.);
	
	fTrueVtxPos[0] = vtx[0];
	fTrueVtxPos[1] = vtx[1];
	fTrueVtxPos[2] = vtx[2];
	fTrueVtxPos[3] = time;
	fTrueVtxPos[4] = 0.;
	
	fTrueVtxPosDouble.push_back(fTrueVtxPos);
}

/************************************************************************************************************************/
// Geometry functions:

double BQFitter::GetDistanceOld(std::vector<double> A, std::vector<double> B) {

	double d1 = (A[0] - B[0]);
	double d2 = (A[0] - B[0]);
	double d3 = (A[0] - B[0]);
	
	return (sqrt(	d1 * d1 +
			d2 * d2 +
			d3 * d3 ) );
}

double GetDistance2(double* A, double* B) {

	double d1 = (A[0] - B[0]);
	double d2 = (A[0] - B[0]);
	double d3 = (A[0] - B[0]);
	
	return (sqrt(	d1 * d1 +
			d2 * d2 +
			d3 * d3 ) );
}

double BQFitter::GroupPMTs(int pmtType, int pmt_number) {
	if(pmtType==0) 			return 0;
	else{
		if(pmt_number<=12) 	return 0;
		else if(pmt_number<=18) return 1;
		else 			return NGROUP_PMT - 1;
	}
}

void BQFitter::LoadPMTInfo() {

	// clear
	fPMT_TubeInMPMT.clear();
	fPMT_Group.clear();
	fPMT_Info.clear();
	/*
	for ( int iPMT=0; iPMT < MAX_PMT; iPMT++ ) {
		fPMT_Info[iPMT][0] 	= 0;
		fPMT_Info[iPMT][1] 	= 0;
		fPMT_Info[iPMT][2]	= 0;
		fPMT_Info[iPMT][3] 	= 0;
		fPMT_Info[iPMT][4]	= 0;
		fPMT_Info[iPMT][5]	= 0;
		
		fPMT_Group[iPMT]	= 0;
		fPMT_TubeInMPMT[iPMT]   = 0;
	}*/
	
	int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
	int iNbr_mPMT = iNbr_Norm;// fWCGeo->GetWCNumPMT(true);
	
	fPMT_Info	.resize(MAX_PMT);
	fPMT_Group	.assign(MAX_PMT,0);
	fPMT_TubeInMPMT	.assign(MAX_PMT,0);
	
	
	WCSimRootPMT wPMT;
	
	// Normal PMTs
	for ( int iPMT=0; iPMT < iNbr_Norm; iPMT++ ) {
	
		wPMT = fWCGeo->GetPMT(iPMT,false);
		std::vector<double> vPMT(6,0);
		
		for ( int j=0; j < 3; j++ ) {
			vPMT[j  ] = wPMT.GetPosition(j);
			vPMT[j+3] = wPMT.GetOrientation(j);
		}
		
		fPMT_TubeInMPMT[iPMT] = 0;
		fPMT_Group     [iPMT] = 0;
		fPMT_Info      [iPMT] = vPMT;
	}
	
	// mPMTs
	for ( int iPMT=0; iPMT < iNbr_mPMT; iPMT++ ) {
	
		wPMT = fWCGeo->GetPMT(iPMT,true);
		std::vector<double> vPMT(6,0);
		
		for ( int j=0; j < 3; j++ ) {
			vPMT[j  ] = wPMT.GetPosition(j);
			vPMT[j+3] = wPMT.GetOrientation(j);
		}
		
		int iNewID = iPMT+mPMT_ID_SHIFT;
		
		fPMT_TubeInMPMT[iNewID] = wPMT.GetmPMT_PMTNo();
		fPMT_Group     [iNewID] = this->GroupPMTs(1,wPMT.GetmPMT_PMTNo());
		fPMT_Info      [iNewID] = vPMT;
	}
}


/************************************************************************************************************************/
// Spline Integral functions:
double BQFitter::SplineIntegral(TSpline3 * s,double start,double end,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*s->Eval(i+stepSize/2);
	}
	return integral;
}
double BQFitter::SplineIntegralAndSubstract(TSpline3 * s0,TSpline3 * s1,double start,double end,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*(s0->Eval(i+stepSize/2)-s1->Eval(i+stepSize/2));
	}
	return integral;
}

double BQFitter::SplineIntegralExpo(TSpline3 * s,double start,double end,double sigma,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*s->Eval(i+stepSize/2)*TMath::Gaus(i+stepSize/2,0,sigma);
	}
	return integral;
}

/************************************************************************************************************************/

double BQFitter::FindDirectionTheta(std::vector<double> vertex,int tubeNumber, int verbose) {

	//clock_t timeStart=clock();
	
	double PMTpos[3];
	double PMTdir[3];
	
	for(int j=0;j<3;j++){
		PMTpos[j] = fPMT_Info[tubeNumber][j  ];
		PMTdir[j] = fPMT_Info[tubeNumber][j+3];
	}

	int pmt_number_in_mpmt = fPMT_TubeInMPMT[tubeNumber];
	double particleRelativePMTpos[3];
	for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vertex[j];
	if(verbose == 3){
		std::cout<<"Original PMT number = "<<pmt_number_in_mpmt<<std::endl;  
		std::cout<<"Direction = "<<PMTdir[0]<<", "<<PMTdir[1]<<", "<<PMTdir[2]<<std::endl;
	}

	double norm=0;
	double ps=0;
	for(int j =0;j<3;j++){
		ps   += (-particleRelativePMTpos[j])*PMTdir[j];
		norm += (particleRelativePMTpos[j] * particleRelativePMTpos[j]);
	}
	double cosTheta=ps/TMath::Sqrt(norm);
	double Theta=TMath::ACos(cosTheta)*180/TMath::Pi();
  
	if(verbose == 3){
		std::cout<<"PMT #"<<pmt_number_in_mpmt<< ", Theta = "<< Theta << std::endl;
	}
  
	//clock_t timeEnd=clock();
	//if(verbose == 3) 
	//	std::cout<<"FindDirectionTheta: Time for the loop = "<<timeEnd-timeStart<<std::endl;
	return Theta;
  
}

void BQFitter::MakeEventInfo(double lowerLimit, double upperLimit) {

	if ( fLastLowerLimit == lowerLimit && fLastUpperLimit == upperLimit ) return;

	fLastLowerLimit = lowerLimit;
	fLastUpperLimit = upperLimit;
	
	for(int i=0; i < NPMT_CONFIGURATION; i++) {
		fEventInfo[i].hits 		= 0;
		fEventInfo[i].SignaloverNoise 	= 0.;
		fEventInfo[i].NoiseIntegral   	= 0.;
		fEventInfo[i].SignalIntegral 	= 0.;
	}
	
	int iHitTotal = fHitInfo.size();
	
	for(int iHit=0; iHit < iHitTotal; iHit++ ) {
		int iPMT = fHitInfo[iHit].PMT;
		int iType = GetPMTType(iPMT);
		fEventInfo[iType].hits += 1;
	}
	
	
	for(int i=0; i < NPMT_CONFIGURATION; i++) {
		double signalDR, signalPE;
		signalDR=fTimeWindowSizeFull*fDarkRate_ns[i];//Over the whole time window
		signalPE=std::max(fEventInfo[i].hits - signalDR,0.);//Over the whole time window.
		
		double signalPETime=signalPE;//I assume that all signal is here.
		double signalDRTime=(upperLimit-lowerLimit)*signalDR/fTimeWindowSizeFull;
		
		fEventInfo[i].SignaloverNoise = signalPETime / signalDRTime;//Over the 
		fEventInfo[i].NoiseIntegral   = this->SplineIntegral(fSplineTimePDFDarkRate[i],lowerLimit,upperLimit);
		fEventInfo[i].SignalIntegral  = this->SplineIntegralAndSubstract(fSplineTimePDFQueue[i],fSplineTimePDFDarkRate[i],lowerLimit,upperLimit);
		
		if(VERBOSE==3) std::cout<<"nhits="<<fEventInfo[i].hits<<", DR average="<<signalDR
					<<", signal over noise="<<fEventInfo[i].SignaloverNoise
					<<", signal integral="  <<fEventInfo[i].SignalIntegral
					<<", DR integral="	<<fEventInfo[i].NoiseIntegral
					<<",in integral="	<<fEventInfo[i].SignalIntegral/fEventInfo[i].NoiseIntegral<<std::endl;
	
		
		//if(VERBOSE==2){
		//	std::cout<<"Lower limit = "<<lowerLimit<<", proba = "<<fSplineTimePDFQueue[i]->Eval(lowerLimit)<<std::endl;
		//	std::cout<<"Upper limit = "<<upperLimit<<", proba = "<<fSplineTimePDFQueue[i]->Eval(upperLimit)<<std::endl;
		//}
	
	}
	

}

void BQFitter::MakePositionList() {
	// Make position list for a given step size

	double dStep = fSearchVtxStep;
	fPositionList.clear();
		
	for(double time = -50; time < 50; time+=(dStep/fLightSpeed)){
		for(double radius = 0; radius <= fTankRadius; radius+=dStep){
		
			double perimeter = 2*TMath::Pi()*radius;
			int numberOfStepsOnCircle = floor(perimeter/dStep);
			if(numberOfStepsOnCircle == 0) numberOfStepsOnCircle = 1;
			double angleStepOnCircle = 2*TMath::Pi()/numberOfStepsOnCircle;
      
			for(double angle=0; angle<=2*TMath::Pi(); angle+=angleStepOnCircle){
			
				for(double height=-fTankHalfHeight;height <= fTankHalfHeight;height+=dStep){
				
					std::vector<double> tVtxPosition(4,0.);//In centimeters
					tVtxPosition[0] = radius*TMath::Cos(angle);
					tVtxPosition[1] = radius*TMath::Sin(angle);
					tVtxPosition[2] = height;
					tVtxPosition[3] = time;
					
					fPositionList.push_back(tVtxPosition);
					
				}
			}
		}
	}
	
	
}

/************************************************************************************************************************/

double BQFitter::FindNLL_Likelihood(std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality){
	
	//TStopwatch timer;
	//timer.Reset();
	//timer.Start();
	
	double NLL = 0;
	
	//double PDFnormalization=1;//this->SplineIntegral(fSplineTimePDFQueue[pmtType],lowerLimit,upperLimit);
	//if(PDFnormalization_fullTimeWindow!=0) PDFnormalization/=PDFnormalization_fullTimeWindow;
	//We know the DR integral over a period of time. But the proba to have DR in one bin depends on the binning of the Spline, since it depends on the bin size.
	//Shouldn't we: 1. Extract the graph from the spline3
	/* TH1D* h = (TH1D*) ->GetHistogram();
	double signalPEinTime=;
	double DRinTime=fDarkRate_ns[pmtType]*(upperLimit-lowerLimit);
	//We should adapt signal / DR the same way as in time window total
	double signalPE = std::max(nhits - fTimeWindowSizeFull*fDarkRate_ns[pmtType],0.);
	double signalIncrease=signalPE/signalPEInTime;//We will scale the signal to have the same integral as number of hits.
	//In that situation, it is simple to deduce the DR value: in a time window 
	double DRtotal = fTimeWindowSizeFull*fDarkRate_ns[pmtType];
	*/
	//double DRtotal = fDarkRate_ns[pmtType]
	//double DR=;//To rescale the PDF, we should know the relative integral of signal vs DR
	
	for(int ihit = 0; ihit < nhits; ihit++){
		int iPMT = fHitInfo[ihit].PMT;
		double hitTime = fHitInfo[ihit].T;
		
		int pmtType = GetPMTType(iPMT);
		double distance = GetDistance(fPMT_Info[iPMT],vertexPosition); 
		
		double tof = distance / fLightSpeed;
		double residual = hitTime - tof - vertexPosition[3];
		double proba = 0;
#ifdef KILLHALF
		if(pmtType==1){
			double t=fRand->Uniform(0,1);
			if(t>0.5) continue;
		}
#endif
		
		bool condition;
		//if(pmtType == 0 || (pmtType == 1 && fLimit_mPMT) ) condition = residual > lowerLimit && residual < upperLimit;
		//else condition = true;
		condition = residual > lowerLimit && residual < upperLimit;
		if(condition){
			//if(pmtType==1) continue;
			proba = fSplineTimePDFQueue[pmtType]->Eval(residual);
			//else proba = 1;//fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
			if(verbose==3){
				std::cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<std::endl;
				std::cout<<"Residual="<<residual<<", proba="<<proba<<", pmt type="<<pmtType<<std::endl;
			}
#endif
			
			if(scaleDR){
				std::cout << "scaleDR with Lower Limit: " << lowerLimit << " Upper Limit: " << upperLimit << std::endl;
				this->MakeEventInfo(lowerLimit,upperLimit);
					
				//First, substract the DR from the PDF to keep only the signal:
				double DR=fSplineTimePDFDarkRate[pmtType]->Eval(residual);
				proba-=DR;
				//Then, scale signal so that signal integral / DR integral = signalOverNoise
				//To do so, we should scale signal so that integral = DR integral * signalOverNoise
				//And of course, to rescale, we should divide by the signal integral
				double Factor = fEventInfo[pmtType].NoiseIntegral * fEventInfo[pmtType].SignaloverNoise / fEventInfo[pmtType].SignalIntegral;
				proba*=Factor;
				//And add again the DR
				proba+=DR;
#ifdef VERBOSE_NLL
				if(verbose==3) std::cout<<"proba after scaling="<<proba<<std::endl;
#endif
			}
			
		}
		else if(killEdges) {
			if ( residual >= upperLimit ){
				//continue;
				proba=fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
				//if(pmtType==1) proba/=1;
#ifdef VERBOSE_NLL
				if(verbose==3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Upper limit = "<<upperLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else if( residual <= lowerLimit ){
				//continue;
				proba=fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
				//if(pmtType==1) proba/=1;
#ifdef VERBOSE_NLL
				if(verbose==3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Lower limit = "<<lowerLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
		}
		else{
			proba=0;
			//cout<<"GUAAAAAAAAAAAAAAAAAAAAAAAH, kill edges="<<killEdges<<", residual="<<residual<<endl;
		}
				
		if(proba<0){
			std::cout<<"Error in PDF"<<std::endl;
			if(proba>-1e-1) proba=0;//Since spline sometimes slightly goes below 0 due to interpolation
			else {
				std::cout <<"Error in "<< residual << "ns where proba = " << proba << std::endl;
				return 0;
			}
		}
		
		if(proba==0) proba=1e-20;
		NLL += -TMath::Log(proba);
	}
	
	//std::cout << " FindNLL_Likelihood loop: " << timer.RealTime() << " " << nhits << std::endl;
	if(directionality!=0){
		//double NLLdirBaye=0;
		double NLLdir=0;
		
		//NLLdirBaye=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,lowerLimit,upperLimit);
		NLLdir=this->FindNLLDirectionality(vertexPosition, nhits, VERBOSE,lowerLimit,upperLimit);
		
#ifdef VERBOSE_NLL
		if(verbose == 2) std::cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdirBaye)<<", dir 2 = "<<NLLdir<<std::endl;
#endif
		if(directionality == 1) NLL+=NLLdir;
		else if(directionality == 2) NLL=NLLdir;
	}
	
	//timer.Stop();
	
	//std::cout << " FindNLL_Likelihood full: " << timer.RealTime() << std::endl;

	return NLL;
}


double BQFitter::FindNLL_NoLikelihood(std::vector<double> vertexPosition, int nhits, double /*lowerLimit*/, double /*upperLimit*/, bool /*killEdges*/, bool /*scaleDR*/, int directionality){
	double NLL = 0;
	
	//double PDFnormalization=1;//this->SplineIntegral(fSplineTimePDFQueue[pmtType],lowerLimit,upperLimit);
	//if(PDFnormalization_fullTimeWindow!=0) PDFnormalization/=PDFnormalization_fullTimeWindow;
	//We know the DR integral over a period of time. But the proba to have DR in one bin depends on the binning of the Spline, since it depends on the bin size.
	//Shouldn't we: 1. Extract the graph from the spline3
	/* TH1D* h = (TH1D*) ->GetHistogram();
	double signalPEinTime=;
	double DRinTime=fDarkRate_ns[pmtType]*(upperLimit-lowerLimit);
	//We should adapt signal / DR the same way as in time window total
	double signalPE = std::max(nhits - fTimeWindowSizeFull*fDarkRate_ns[pmtType],0.);
	double signalIncrease=signalPE/signalPEInTime;//We will scale the signal to have the same integral as number of hits.
	//In that situation, it is simple to deduce the DR value: in a time window 
	double DRtotal = fTimeWindowSizeFull*fDarkRate_ns[pmtType];
	*/
	//double DRtotal = fDarkRate_ns[pmtType]
	//double DR=;//To rescale the PDF, we should know the relative integral of signal vs DR
			
	//std::cout << " Find NLL " << nhits << std::endl;
	//std::cout << " First Hit " << fHitInfo[0].PMT << " " << vertexPosition.size() << std::endl;	
	for(int ihit = 0; ihit < nhits; ihit++){
	
		//std::cout << " NLL Hit " << ihit << " " << fHitInfo[ihit].PMT << std::endl;
		int iPMT = fHitInfo[ihit].PMT;		
	
		double hitTime = fHitInfo[ihit].T;
		
		int pmtType = GetPMTType(iPMT);
		double distance = GetDistance(fPMT_Info[iPMT],vertexPosition);
		
		double tof = distance / fLightSpeed;
		double residual = hitTime - tof - vertexPosition[3];
#ifdef KILLHALF
		if(pmtType==1){
			double t=fRand->Uniform(0,1);
			if(t>0.5) continue;
		}
#endif

		bool bCondition = (residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive) || (pmtType == 1 && !fLimit_mPMT);

      		//std::cout << " HIT " << ihit << " Has " <<   bCondition << " " << pmtType << " " << fLimit_mPMT << " " << residual << " ( " << hitTime << " - " << tof << " - " << vertexPosition[3] << " ) "  << distance << " / " << fLightSpeed << " " << fPMT_Info[iPMT][0] << " "<< fPMT_Info[iPMT][1] << " "<< fPMT_Info[iPMT][2] << " " << iPMT << std::endl;
      
		if( bCondition ){
			NLL++;//= fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
			if(verbose==3){
				std::cout<<"Residual="<<residual<<", pmt type="<<pmtType<<std::endl;
				std::cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<std::endl;
				std::cout<<"residual="<<residual<<std::endl;
			}
#endif
		}
	}
	
	//std::cout << " NLL Final" << NLL << std::endl;
	
	//std::cout<<"NLL="<<NLL<<std::endl;
	//if(NLL!=0) NLL=1/NLL;
	
	if(NLL!=0) NLL=-TMath::Log(NLL);
	else NLL=1e15;
	
	//if(verbose==2) cout<<"NLL="<<NLL<<endl;
	if(directionality!=0){
		//double NLLdirBaye=0;
		double NLLdir=0;
		
		//NLLdirBaye=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,fHitTimeLimitsNegative,fHitTimeLimitsPositive);
		NLLdir=this->FindNLLDirectionality(vertexPosition, nhits, VERBOSE,fHitTimeLimitsNegative,fHitTimeLimitsPositive);
#ifdef VERBOSE_NLL
			
		if(verbose == 2) std::cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdirBaye)<<", dir = "<<NLLdir<<std::endl;
#endif
		if(directionality == 1) NLL+=NLLdir;
		else if(directionality == 2) NLL=NLLdir;
	}
	
	//std::cout<<"NLL="<<NLL<<std::endl;

	return NLL;
}


double BQFitter::FindNLL(std::vector<double> vertexPosition, int nhits, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality){
	double NLL = 0;
	
	//double PDFnormalization=1;//this->SplineIntegral(fSplineTimePDFQueue[pmtType],lowerLimit,upperLimit);
	//if(PDFnormalization_fullTimeWindow!=0) PDFnormalization/=PDFnormalization_fullTimeWindow;
	//We know the DR integral over a period of time. But the proba to have DR in one bin depends on the binning of the Spline, since it depends on the bin size.
	//Shouldn't we: 1. Extract the graph from the spline3
	/* TH1D* h = (TH1D*) ->GetHistogram();
	double signalPEinTime=;
	double DRinTime=fDarkRate_ns[pmtType]*(upperLimit-lowerLimit);
	//We should adapt signal / DR the same way as in time window total
	double signalPE = std::max(nhits - fTimeWindowSizeFull*fDarkRate_ns[pmtType],0.);
	double signalIncrease=signalPE/signalPEInTime;//We will scale the signal to have the same integral as number of hits.
	//In that situation, it is simple to deduce the DR value: in a time window 
	double DRtotal = fTimeWindowSizeFull*fDarkRate_ns[pmtType];
	*/
	//double DRtotal = fDarkRate_ns[pmtType]
	//double DR=;//To rescale the PDF, we should know the relative integral of signal vs DR
		
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	
	for(int ihit = 0; ihit < nhits; ihit++){
		int iPMT = fHitInfo[ihit].PMT;
		double hitTime = fHitInfo[ihit].T;
		
		int pmtType = GetPMTType(iPMT);
		double distance = GetDistance(fPMT_Info[iPMT],vertexPosition);
		
		double tof = distance / fLightSpeed;
		double residual = hitTime - tof - vertexPosition[3];
		double proba;
#ifdef KILLHALF
		if(pmtType==1){
			double t=fRand->Uniform(0,1);
			if(t>0.5) continue;
		}
#endif

		if(likelihood){
			bool condition;
			//if(pmtType == 0 || (pmtType == 1 && fLimit_mPMT) ) condition = residual > lowerLimit && residual < upperLimit;
			//else condition = true;
			condition = residual > lowerLimit && residual < upperLimit;
			if(condition){
				//if(pmtType==1) continue;
				proba = fSplineTimePDFQueue[pmtType]->Eval(residual);
				//else proba = 1;//fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if(verbose==3){
					std::cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<std::endl;
					std::cout<<"Residual="<<residual<<", proba="<<proba<<", pmt type="<<pmtType<<std::endl;
				}
#endif
				if(scaleDR){
					std::cout << "scaleDR with Lower Limit: " << lowerLimit << " Upper Limit: " << upperLimit << std::endl;
					this->MakeEventInfo(lowerLimit,upperLimit);
					
					//First, substract the DR from the PDF to keep only the signal:
					double DR=fSplineTimePDFDarkRate[pmtType]->Eval(residual);
					proba-=DR;
					//Then, scale signal so that signal integral / DR integral = signalOverNoise
					//To do so, we should scale signal so that integral = DR integral * signalOverNoise
					//And of course, to rescale, we should divide by the signal integral
					double Factor = fEventInfo[pmtType].NoiseIntegral * fEventInfo[pmtType].SignaloverNoise / fEventInfo[pmtType].SignalIntegral;
					proba*=Factor;
					//And add again the DR
					proba+=DR;
					if(verbose==3) std::cout<<"proba after scaling="<<proba<<std::endl;
				}
			}
			else if(killEdges && residual >= upperLimit){
				//continue;
				proba=fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
				//if(pmtType==1) proba/=1;
#ifdef VERBOSE_NLL
				if(verbose==3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Upper limit = "<<upperLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else if(killEdges && residual <= lowerLimit){
				//continue;
				proba=fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
				//if(pmtType==1) proba/=1;
#ifdef VERBOSE_NLL
				if(verbose==3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Lower limit = "<<lowerLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else{
				proba=0;
				//cout<<"GUAAAAAAAAAAAAAAAAAAAAAAAH, kill edges="<<killEdges<<", residual="<<residual<<endl;
			}
		}
		else{
			bool condition;
			if(pmtType == 0 || (pmtType == 1 && fLimit_mPMT) ) condition = residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive;
			else condition = true;
			if(condition){
				NLL++;//= fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if(verbose==3){
					std::cout<<"Residual="<<residual<<", pmt type="<<pmtType<<std::endl;
					std::cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<std::endl;
					std::cout<<"residual="<<residual<<", proba="<<proba<<std::endl;
				}
#endif
			}
			else if(killEdges) continue;
			else proba=0;
		}

		if(likelihood){
			if(proba<0){
				std::cout<<"Error in PDF"<<std::endl;
				if(proba>-1e-1) proba=0;//Since spline sometimes slightly goes below 0 due to interpolation
				else {
					std::cout <<"Error in "<< residual << "ns where proba = " << proba << std::endl;
					return 0;
				}
			}
			if(proba==0) proba=1e-20;
			NLL += -TMath::Log(proba);
		}
	}
	std::cout<<"FindNLL: Time it took for the loop = "<< timer.RealTime() << std::endl;
	//cout<<"NLL="<<NLL<<endl;
	if(!likelihood){
		//if(NLL!=0) NLL=1/NLL;
		if(NLL!=0) NLL=-TMath::Log(NLL);
		else NLL=1e15;
	}
	//if(verbose==2) cout<<"NLL="<<NLL<<endl;
	if(directionality!=0){
		double NLLdir=0;
		double NLLdir2=0;
		if(likelihood){
			//NLLdir=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,lowerLimit,upperLimit);
			NLLdir2=this->FindNLLDirectionality(vertexPosition, nhits, verbose,lowerLimit,upperLimit);
		}
		else{
			//NLLdir=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,fHitTimeLimitsNegative,fHitTimeLimitsPositive);
			NLLdir2=this->FindNLLDirectionality(vertexPosition, nhits, verbose,fHitTimeLimitsNegative,fHitTimeLimitsPositive);
		}
		if(verbose == 2) std::cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdir)<<", dir 2 = "<<NLLdir2<<std::endl;
		if(directionality == 1) NLL+=NLLdir2;
		else if(directionality == 2) NLL=NLLdir2;
	}
	timer.Stop();
	std::cout<<"FindNLL: Total time = "<< timer.RealTime() << std::endl;
	timer.Reset();
	

	return NLL;
}

double BQFitter::FindNLLDirectionality(std::vector<double> vVtxPos, int nhits, int verbose, double lowerLimit, double upperLimit) {

	double NLL = 0;
	//double PDFnormalization=1;//this->SplineIntegral(fSplineTimePDFQueue[pmtType],lowerLimit,upperLimit);
	std::vector<double> vDirection(2,0.);//Return phi and theta.

	for(int ihit = 0; ihit < nhits; ihit++){
		int iPMT = fHitInfo[ihit].PMT;
		double hitTime = fHitInfo[ihit].T;
		int pmtType = GetPMTType(iPMT);
		int tubeNumber = iPMT;
		//int pmt_number_in_mpmt=fPMT_TubeInMPMT[iPMT];

		if(pmtType==0) continue;
		    
		double distance = GetDistance(fPMT_Info[iPMT],vVtxPos);
		double tof = distance / fLightSpeed;
		double residual = hitTime - tof - vVtxPos[3];
		bool condition;
		if(pmtType == 0 || (pmtType == 1 && fLimit_mPMT) ) condition = residual > lowerLimit && residual < upperLimit;
		else condition = true;
		if(condition){
			//if(1/*residual > lowerLimit && residual < upperLimit*/){//Do not use it as fluctuating the number of hits used in one event does allow a fair comparison with this likelihood. If we want to restrict the likelihood, use the "DirectionalNLLBayes".
			//vectorVertexPMT(vPos,tubeNumber,pmtType,vDirection,pmt_number_in_mpmt,verbose);
			//proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(vDirection[1]));
			double theta = this->FindDirectionTheta(vVtxPos,tubeNumber,verbose);
			int pmtGroup = fPMT_Group[tubeNumber];
			double proba;
			proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(theta));

      
			if(proba==0){
				/*
				double min=0;
				for(int ibinx=hPMTDirectionality_1D[pmtType][pmtGroup]->GetNbinsX();ibinx>=1;ibinx--){//Search the last bin of the distri non-zero and use its value
					if(hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx) !=0 ){
						min = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx);
						break;
					}
				}
				//proba=min;
				*/
				proba=1e-20;
				if(verbose) std::cout<<"We are at proba = 0, theta = "<<theta<<", proba used = "<<proba<<std::endl;
			}
			NLL += -TMath::Log(proba);
      
			if(verbose==3){
				//std::cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << vDirection[1] << ", proba="<<proba<<std::endl;
				std::cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << theta << ", proba="<<proba<<std::endl;
			}
		}
	}
	if(verbose==2) std::cout<<"NLL directionnel="<<NLL<<std::endl;
	return NLL;
}
/************************************************************************************************************************/

//Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
//The step between each grid points is given by stepSize (in centimeters).
//We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
//The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower). 
//The number of output vertices can be chosen using "tolerance".
std::vector< std::vector<double> > BQFitter::SearchVertex(int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	
	std::vector<double> 			tBestReconstructedVertexPosition;
	
	//double ** reconstructedVertexPosition = new double[4];
		
	std::vector<struct FitPosition>  tVtxContainer;
	
	clock_t timeStart=clock();
		
	if ( fPositionList.size() == 0 ) std::cout << " Position list not filled " << std::endl;
	
	for ( unsigned int iPos = 0; iPos < fPositionList.size(); iPos++ ) {
		
  		struct FitPosition hPos;
  		hPos.Vtx = fPositionList[iPos];
  		hPos.NLL = 0;
  							
#ifdef CHECK_TO_TRUE_VTX
		double distanceTrue= GetDistance(hPos.Vtx,fTrueVtxPos);
#ifdef VERBOSE_VTX
		if(verbose == 2) std::cout<<"time distance to true = "<< hPos.Vtx[3] - fTrueVtxPos[3] <<"ns, distance to true = " << distanceTrue << "m, "<<std::endl;
#endif
#endif
		if ( likelihood ) 	hPos.NLL = this->FindNLL_Likelihood  (hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		else 			hPos.NLL = this->FindNLL_NoLikelihood(hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		
		hPos.NLL += fRand->Uniform(0,1e-5);//Here only not to have exactly the same value for 2NLL, and be removed from map
	  
	  	tVtxContainer.push_back(hPos);
	  	
		//std::cout << "Test Pos: " << hPos.Vtx[3] << " " << hPos.Vtx[0] << " " << hPos.Vtx[1] << " " << hPos.Vtx[2] << " NLL = " << hPos.NLL << " " << likelihood << std::endl;
					
#ifdef VERBOSE_VTX
		//if(verbose && distanceToTrue < 2) std::cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<std::endl;
#endif
	}
	
	std::sort(tVtxContainer.begin(),tVtxContainer.end(),SortingNLL());
	while((int) tVtxContainer.size() > tolerance) {		
#ifdef VERBOSE_VTX
		std::cout << " Remove end: " << tVtxContainer.back().NLL << " (first " << tVtxContainer[0].NLL << ") " << tVtxContainer.size() << std::endl;
#endif
		tVtxContainer.pop_back();
	}
	
	clock_t timeEndLoop=clock();
	
	std::vector< std::vector<double> > tReconstructedVertexPosition;
	for ( int iPos = 0; iPos < tolerance; iPos++ ) {
		std::vector<double> vPos(5,0.);
		struct FitPosition hPos = tVtxContainer[iPos];
		
		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;
		
		tReconstructedVertexPosition.push_back(vPos);
		
#ifdef VERBOSE_VTX
		double distanceToTrue = GetDistance(vPos,fTrueVtxPos);
			
		std::cout<<"NLL = " << hPos.NLL << ", time distance to true = " << vPos[3]-fTrueVtxPos[3] << "ns, distance to true = " << distanceToTrue <<"m, "<<std::endl;
		std::cout<<"NLL = " << hPos.NLL << ", vertex = " << vPos[0]<<", "<<vPos[1]<<", "<<vPos[2]<<", "<<vPos[3]<<std::endl;
#endif
		std::cout<<"NLL = " << hPos.NLL << ", vertex = " << vPos[0]<<", "<<vPos[1]<<", "<<vPos[2]<<", "<<vPos[3]<<std::endl;
	}
	
#ifdef CHECK_TO_TRUE_VTX
	tBestReconstructedVertexPosition = tReconstructedVertexPosition[0]; 
	
	double distanceToTrue = GetDistance(tBestReconstructedVertexPosition,fTrueVtxPos);
#ifdef VERBOSE_VTX
	std::cout << "Candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << tBestReconstructedVertexPosition[3] - fTrueVtxPos[3] << "ns, with a min NLL = " <<minNLL << std::endl;
	std::cout << "Time of the true event = " << fTrueVtxPos[3] << std::endl;
#endif
#endif
	clock_t timeEnd=clock();
	
	//if(VERBOSE) 
		std::cout<<"SearchVertex: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<< " " << fPositionList.size() << std::endl;
	
	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
}

std::vector< std::vector<double> > BQFitter::SearchVertex_Main(	
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
		
		
	std::vector<std::thread> lThreadList;
	std::vector< std::vector<double> > lOutputFinal;
	
	int iCand_Step = (int) fPositionList.size() / N_THREAD;
	
	if ( iCand_Step < 1 ) iCand_Step = 1;
	
	for ( int iStart = 0; iStart < (int) fPositionList.size(); iStart += iCand_Step ) {
	
		
		std::thread tThrd( SearchVertex_CallThread, iStart,iCand_Step, 
							nhits, tolerance, likelihood, lowerLimit, upperLimit, directionality);
							
		lThreadList.push_back(std::move(tThrd));
		
	}
	
	
	for ( unsigned int iIdx=0; iIdx < lThreadList.size(); iIdx++ ) {
		// Wait thread
		
		if ( lThreadList[iIdx].joinable() )
			lThreadList[iIdx].join();
	}
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	lOutputFinal = fThreadOutput;
	fThreadOutput.clear();
	mtx.unlock();
		
	std::sort(lOutputFinal.begin(), lOutputFinal.end(), SortOutputVector); 	

	while((int) lOutputFinal.size() > tolerance) {	
		lOutputFinal.pop_back();
	}
	
	//std::cout << " SearchVtx main " << tolerance << " " << lOutputFinal.size() << std::endl;
	return lOutputFinal;
}

//Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
//The step between each grid points is given by stepSize (in centimeters).
//We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
//The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower). 
//The number of output vertices can be chosen using "tolerance".
void BQFitter::SearchVertex_thread(
			int iStart, int iIte,
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	
	std::vector<struct FitPosition>  tVtxContainer;
			
	if ( fPositionList.size() == 0 ) std::cout << " Position list not filled " << std::endl;
	
	int iEnd = (iIte+iStart);
	if ( iEnd > (int) fPositionList.size() ) iEnd = fPositionList.size();
	
	for ( int iPos = iStart; iPos < iEnd; iPos++ ) {
				
  		struct FitPosition hPos;
  		hPos.Vtx = fPositionList[iPos];
  		hPos.NLL = 0;
		//std::cout << " Thread2 " << iPos << " " << iPos-iStart << " " << iEnd <<  std::endl;	
  		
		if ( likelihood ) 	hPos.NLL = this->FindNLL_Likelihood  (hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		else 			hPos.NLL = this->FindNLL_NoLikelihood(hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		
		//std::cout << " Thread3 " << iPos << " " << iPos-iStart << " " << iEnd <<  std::endl;	
		hPos.NLL += fRand->Uniform(0,1e-5);//Here only not to have exactly the same value for 2NLL, and be removed from map
	  
	  	tVtxContainer.push_back(hPos);
	}
	
	std::sort(tVtxContainer.begin(),tVtxContainer.end(),SortingNLL());
	while((int) tVtxContainer.size() > tolerance) {		
		tVtxContainer.pop_back();
	}
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	//std::vector< std::vector<double> > tReconstructedVertexPosition;
	for ( unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++ ) {
		std::vector<double> vPos(5,0.);
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


//Fine grid search of the vertex. It will search in a box this time, and not a cylinder as searchVertex is doing. It is useful when we have one or several first candidate vertices. Therefore, this code is generally run after serachVertex.
//The box is centered around the initial vertex provided.
//The box size is 2 x limits in each direction, centered on the initial vertex provided.
//nCandidate provide the size of the list of inital vertices, while tolerance provide the size of the output list.
//Other informations are the same as searchVertex.
std::vector< std::vector<double> > BQFitter::SearchVertexFine(std::vector< std::vector<double> > initialVertex,double * limits, double stepSize,int nhits,int nCandidates,int tolerance,int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
	
	if(verbose) std::cout<<"Start fine search, use directionality? "<<directionality<<std::endl;
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	//double * reconstructedVertexPosition = new double[4];
	std::vector<double> tBestReconstructedVertexPosition;
	double minNLL=1e15;
	clock_t timeStart=clock();

	std::map< double, std::vector<double> > vertexContainer;
  
	for(int icand=0;icand<nCandidates;icand++){
		if(verbose) std::cout<<"Candidate #"<<icand<<std::endl;
		for(double time = initialVertex[icand][VTX_T] - limits[0]/fLightSpeed; time <= initialVertex[icand][VTX_T]+limits[0]/fLightSpeed; time+=(stepSize/fLightSpeed)){
			if(verbose == 2) std::cout << "Time = " << time << "ns" << std::endl;
			for(double x = initialVertex[icand][VTX_X] - limits[1]; x <= initialVertex[icand][VTX_X] + limits[1];x+=stepSize){
				if(verbose == 2) std::cout << "x = " << x << "m" << std::endl;
				for(double y = initialVertex[icand][VTX_Y] - limits[2]; y <= initialVertex[icand][VTX_Y] + limits[2];y+=stepSize){
					if(verbose == 2) std::cout << "y = " << y << "m" << std::endl;
					for(double z = initialVertex[icand][VTX_Z] - limits[3]; z <= initialVertex[icand][VTX_Z] + limits[3];z+=stepSize){
						if(verbose == 2) std::cout << "z = " << z << "m" << std::endl;

						std::vector<double> tVtxPosition(4,0.);//In centimeters
						tVtxPosition[0] = x;//radius*TMath::Cos(angle);
						tVtxPosition[1] = y;//radius*TMath::Sin(angle);
						tVtxPosition[2] = z;//height;
						tVtxPosition[3] = time;

						double distanceTrue = GetDistance(tVtxPosition,fTrueVtxPos);
						if(verbose == 2) std::cout<<"time distance to true = "<<tVtxPosition[3]-fTrueVtxPos[3]<<"ns, distance to true = "<<distanceTrue<<"m, "<<std::endl;

						double NLL=0;
						if(average){
							if(verbose==2) std::cout << std::endl << "Vertex timing = "<<tVtxPosition[3]-fTrueVtxPos[3]<<", (x,y,z)=("<<tVtxPosition[0]-fTrueVtxPos[0]<<","<<tVtxPosition[1]-fTrueVtxPos[1]<<","<<tVtxPosition[2]-fTrueVtxPos[2]<<")"<<std::endl;
							//double vNLL[fAveraging];
							for(int irand=0;irand<fAveraging;irand++){
								//Determine the position around the vertex
								double radius = 1.5*fRand->Uniform(0,stepSize);
								double phi = fRand->Uniform(0,2*TMath::Pi());
								double theta = fRand->Uniform(0,TMath::Pi());
								double timeRand = 1.5*fRand->Uniform(-stepSize/fLightSpeed,stepSize/fLightSpeed);
								std::vector<double> randVertexPosition(4,0.);//In centimeters
								randVertexPosition[0] = tVtxPosition[0] + radius*TMath::Sin(theta)*TMath::Cos(phi);
								randVertexPosition[1] = tVtxPosition[1] + radius*TMath::Sin(theta)*TMath::Sin(phi);
								randVertexPosition[2] = tVtxPosition[2] + radius*TMath::Cos(theta);
								randVertexPosition[3] = tVtxPosition[3] + timeRand;//rand->Uniform(-stepSize/fLightSpeed/2,stepSize/fLightSpeed/2);
								//Evaluate its NLL
								double nll=this->FindNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);
								//double nll=this->FindNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit);
								//double nll=findRMS(iPMTConfiguration,randVertexPosition,fHitInfo,nhits,verbose,lowerLimit,upperLimit);
								if(verbose==2) std::cout << "Randomized Vertex timing = "<<randVertexPosition[3]-fTrueVtxPos[3]<<", (x,y,z)=("<<randVertexPosition[0]-fTrueVtxPos[0]<<","<<randVertexPosition[1]-fTrueVtxPos[1]<<","<<randVertexPosition[2]-fTrueVtxPos[2]<<"), nll="<<nll<<std::endl;
								NLL+=nll;
								//vNLL[irand]=nll;
							}
	      
							if(fAveraging!=0) NLL/=fAveraging;
							//double rmsNLL=0;
							//for(int irand=0;irand<fAveraging;irand++){
							//rmsNLL+=pow(vNLL[irand]-NLL,2);
							//}
							//if(fAveraging!=0) rmsNLL/=fAveraging;
							//rmsNLL=TMath::Sqrt(rmsNLL);
							//NLL=rmsNLL;
						}
						else NLL=this->FindNLL(tVtxPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);
						//else NLL=findChi2BinnedExpo(iPMTConfiguration,tVtxPosition,nhits,1,100);
						//else NLL=findRMS(iPMTConfiguration,tVtxPosition,fHitInfo,nhits,verbose,lowerLimit,upperLimit);

	    
						std::vector<double> candidateReconstructedVertexPosition(4,0.);
						for(int i=0;i<4;i++) candidateReconstructedVertexPosition[i] = tVtxPosition[i];
						//for(int i=0;i<4;i++) candidateReconstructedVertexPosition.push_back(vertexPosition[i]);
						//reconstructedVertexPositionArray.push_back(candidateReconstructedVertexPosition)//;
						//std::map< double, std::vector<double> >::iterator itMap = vertexContainer.end()-1;
	    
						vertexContainer[NLL] = candidateReconstructedVertexPosition;
						if( (int) vertexContainer.size() > tolerance){
							vertexContainer.erase(--vertexContainer.rbegin().base());
						}
						/*
						std::vector<double> candidateReconstructedVertexPosition;//[4];// = new double[4];
						for(int i=0;i<4;i++) candidateReconstructedVertexPosition.push_back(vertexPosition[i]);
						//reconstructedVertexPositionArray.push_back(candidateReconstructedVertexPosition);
						vertexContainer[NLL] = candidateReconstructedVertexPosition;
						*/
						if(NLL<minNLL){
							minNLL=NLL;
							//std::cout << "minNLL = "<<minNLL << std::endl;
							tBestReconstructedVertexPosition = tVtxPosition;
						}/*
						candidateReconstructedVertexPosition.clear();*/
						//if(verbose && distanceToTrue < 2) std::cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<std::endl;
						//if(verbose && y=) std::cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<std::endl;
					}
				}
			}
		}
	}
	clock_t timeEndLoop=clock();
	std::vector< std::vector<double> > tReconstructedVertexPosition;
	for(int count=0;count<tolerance;count++){
		std::vector<double> vTmp(5,0.);
		tReconstructedVertexPosition.push_back(vTmp);
	}
	int counter=0;
	for(std::map< double, std::vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
		//for(std::map< double, std::vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
		if(counter<tolerance){
			//tReconstructedVertexPosition[counter] = new double[4];
			for(int i=0;i<4;i++) tReconstructedVertexPosition[counter][i] = (it->second)[i];
				tReconstructedVertexPosition[counter][4] = it->first;

				if(verbose){
					double distanceToTrue=0;
					for(int i=0;i<3;i++){
						distanceToTrue += pow(fTrueVtxPos[i]-it->second[i],2);
					}
					distanceToTrue = TMath::Sqrt(distanceToTrue);
					std::cout<<"NLL = "<<it->first<<", time distance to true = "<<it->second[3]-fTrueVtxPos[3]<<"ns, distance to true = "<<distanceToTrue<<"m, "<<std::endl;
					//if(verbose) std::cout<<"NLL = "<<it->first<<", vertex = "<<it->second[0]<<", "<<it->second[1]<<", "<<it->second[2]<<", "<<it->second[3]<<std::endl;
				}
			}
    
		counter++;
	}
	double distanceToTrue=0;
	for(int i=0;i<3;i++){
		distanceToTrue += pow(fTrueVtxPos[i]-tBestReconstructedVertexPosition[i],2);
	}
	distanceToTrue = TMath::Sqrt(distanceToTrue);
	if(verbose){
		std::cout << "Fine candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << tBestReconstructedVertexPosition[3] - fTrueVtxPos[3] << "ns, with a min NLL = " <<minNLL << std::endl;
		std::cout << "Time of the true event = " << fTrueVtxPos[3] << std::endl;
	}
	clock_t timeEnd=clock();
	if(verbose) 
		std::cout<<"SearchVertexFine: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<std::endl<<std::endl;
  
	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
  
}


//here is the function where minimization gonna take place. The loop is inside this function
std::vector< std::vector<double> > BQFitter::MinimizeVertex(std::vector< std::vector<double> > initialVertex,double * limits, double stepSize,int nhits,int nCandidates,int tolerance,int verbose,bool /*likelihood*/,bool /*average*/,double lowerLimit, double upperLimit, int directionality){
	if(verbose==2) std::cout<<"Minimizer"<<std::endl;
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	//double * reconstructedVertexPosition = new double[4];
	std::vector<double> tBestReconstructedVertexPosition;
	clock_t timeStart=clock();

	std::vector<struct FitPosition>  tVtxContainer;

	if(verbose==2) std::cout<<"Minimizer"<<std::endl;
	TFitter * minimizer = new TFitter(8);//4=nb de params?
	TMinuit * minuit = minimizer->GetMinuit();
	
	double arglist[20];
	int err=0;
	//arglist[0]=0;
	double p1=verbose-1;
	minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
	minuit->SetErrorDef(1);
	minuit->SetPrintLevel(-1); // quiet mode
	arglist[0]=2;
	minuit->mnexcm("SET STR",arglist,1,err);//set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
	minimizer->SetFCN(MinuitLikelihood);//here is function to minimize

	double Mig[2]={1e6,1e0};//maxcalls and tolerance
	//double Imp2[1]={10000};
	
	minimizer->SetParameter(0,"vertex0",0,stepSize		  ,-limits[0]			,limits[0]		);
	minimizer->SetParameter(1,"vertex1",0,stepSize		  ,-limits[1]			,limits[1]		);
	minimizer->SetParameter(2,"vertex2",0,stepSize		  ,-limits[2]			,limits[2]		);
	minimizer->SetParameter(3,"vertex3",0,stepSize/fLightSpeed,-limits[3]/fLightSpeed	,limits[3]/fLightSpeed	);
	
	minimizer->SetParameter(4,"PMTConfiguration",0,0,0,0); // Useless now
	minimizer->SetParameter(5,"nhits",nhits,nhits,nhits-1,nhits+1);
	minimizer->SetParameter(6,"lowerLimit",lowerLimit,5e-1,-7,-2);
	minimizer->SetParameter(7,"upperLimit",upperLimit,5e-1,2,10);
	//minimizer->SetParameter(8,"scaling factor",1,1,0,1e3);
	//minimizer->SetParameter(8,"signal pe",1,1,0,1e3);
	minimizer->SetParameter(8,"expoSigma",100,5,0,1e3);
	minimizer->SetParameter(9,"directionality",directionality,1,0,2);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);

	for(int icand=0;icand<nCandidates;icand++){
#ifdef VERBOSE_MIN
		if(verbose>=1){
			std::cout<<"Candidate #"<<icand<<std::endl;
			std::cout<<"Initial (x,y,z,t): " 	<< initialVertex[icand][0] 
							<<", " 	<< initialVertex[icand][1]
							<<", " 	<< initialVertex[icand][2]
							<<", " 	<< initialVertex[icand][3] << std::endl;
		}
#endif
								
		minimizer->SetParameter(0,"vertex0",initialVertex[icand][0],stepSize		  ,initialVertex[icand][0]-limits[0]			,initialVertex[icand][0]+limits[0]		);
		minimizer->SetParameter(1,"vertex1",initialVertex[icand][1],stepSize		  ,initialVertex[icand][1]-limits[1]			,initialVertex[icand][1]+limits[1]		);
		minimizer->SetParameter(2,"vertex2",initialVertex[icand][2],stepSize		  ,initialVertex[icand][2]-limits[2]			,initialVertex[icand][2]+limits[2]		);
		minimizer->SetParameter(3,"vertex3",initialVertex[icand][3],stepSize/fLightSpeed  ,initialVertex[icand][3]-limits[3]/fLightSpeed	,initialVertex[icand][3]+limits[3]/fLightSpeed	);
	
		minimizer->SetParameter(4,"PMTConfiguration",0,0,0,0); // Useless now
		minimizer->SetParameter(5,"nhits",nhits,nhits,nhits-1,nhits+1);
		minimizer->SetParameter(6,"lowerLimit",lowerLimit,5e-1,-7,-2);
		minimizer->SetParameter(7,"upperLimit",upperLimit,5e-1,2,10);
		minimizer->SetParameter(8,"expoSigma",100,5,0,1e3);
		minimizer->SetParameter(9,"directionality",directionality,1,0,2);
		//minimizer->SetParameter(8,"scaling factor",1,1,0,1e3);
		//minimizer->SetParameter(8,"signal pe",1,1,0,1e3);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);
      
		//double nll=findChi2BinnedExpo(iPMTConfiguration,initialVertex[icand],nhits,1,100);
		//double nll=this->FindNLL(initialVertex[icand],nhits,likelihood,2,lowerLimit,upperLimit,true,true);
		//double nll2=this->FindNLL(initialVertex[icand],nhits,likelihood,2,-3,4);
		//std::cout<<"NLL="<<nll<<std::endl;
		/*
		TFile * fScan = new TFile("fScan.root","recreate");
		TGraph *gr;
		double arglist1[2];
		for(int ipar=0;ipar<4;ipar++){
			arglist1[0]=ipar+1;//param to scan
			arglist1[1]=10000;//number of points
			//arglist1[2]=-200;//lower limit of the scan
			//arglist1[3]=200;//upper limit
			minimizer->ExecuteCommand("SCAN",arglist1,2);//arglist,1);//minmize log like fitting parameters
			minuit=minimizer->GetMinuit();
			gr = (TGraph*)minuit->GetPlot();
			gr->Write(Form("scan_par%d",ipar));
		}
		fScan->Close();
		*/
		//initialVertex[icand][0] = 6.40055e-01;
		//initialVertex[icand][1] = 3.26636e+03;
		//initialVertex[icand][2] = -1.14015e+03;
		//initialVertex[icand][3] = 8.94242e+02;
		//double nll=findChi2Binned(iPMTConfiguration,initialVertex[icand],nhits,1,lowerLimit,upperLimit);
		//minuit->mnseek();
		//minimizer->ExecuteCommand("SIMPLEX",Mig,2);
		/* minimizer->FixParameter(1);
		minimizer->FixParameter(2);
		minimizer->FixParameter(3);
		minimizer->ExecuteCommand("MIGRAD",Mig,2);
		minimizer->ReleaseParameter(1);
		minimizer->ReleaseParameter(2);
		minimizer->ReleaseParameter(3);
		minimizer->FixParameter(0);
		minimizer->ExecuteCommand("MIGRAD",Mig,2);
		minimizer->ReleaseParameter(0);
		minimizer->ExecuteCommand("MIGRAD",Mig,2);*/
		minimizer->ExecuteCommand("MIGRAD",Mig,2);
		//minimizer->ExecuteCommand("SEEK",Mig,2);
		//minimizer->ExecuteCommand("MINIMIZE",Mig,2);
		//minimizer->ExecuteCommand("MINOS",Mig,2);


		std::vector<double> tCandRecoVtxPos(4,0.);
		
		tCandRecoVtxPos[0] = minimizer->GetParameter(0);
		tCandRecoVtxPos[1] = minimizer->GetParameter(1);
		tCandRecoVtxPos[2] = minimizer->GetParameter(2);
		tCandRecoVtxPos[3] = minimizer->GetParameter(3);
		
		struct FitPosition hPos;
		hPos.Vtx = tCandRecoVtxPos;
		hPos.NLL = 0;
		
		// Get minuit
		minuit=minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar,nparx,istat;
		minuit->mnstat(fmin,fedm,errdef,npar,nparx,istat);
		
		hPos.NLL = fmin;
		
#ifdef VERBOSE_MIN
		if(verbose>=1) std::cout<<"NLL="<<hPos.NLL<<std::endl;
#endif
		tVtxContainer.push_back(hPos);
	}	 
	 
	std::sort(tVtxContainer.begin(),tVtxContainer.end(),SortingNLL());
	while((int) tVtxContainer.size() > tolerance) {	
		tVtxContainer.pop_back();
	}
	
	clock_t timeEndLoop=clock();
	
	std::vector< std::vector<double> > tReconstructedVertexPosition;
	for ( int iPos = 0; iPos < tolerance; iPos++ ) {
		std::vector<double> vPos(5,0.);
		struct FitPosition hPos = tVtxContainer[iPos];
		
		vPos[0] = hPos.Vtx[0];
		vPos[1] = hPos.Vtx[1];
		vPos[2] = hPos.Vtx[2];
		vPos[3] = hPos.Vtx[3];
		vPos[4] = hPos.NLL;
		
		tReconstructedVertexPosition.push_back(vPos);
		
#ifdef CHECK_TO_TRUE_VTX
#ifdef VERBOSE_MIN
		double distanceToTrue = GetDistance(vPos,fTrueVtxPos);
			
		std::cout<<"NLL = " << hPos.NLL << ", time distance to true = " << vPos[3]-fTrueVtxPos[3] << "ns, distance to true = " << distanceToTrue <<"m, "<<std::endl;
		std::cout<<"NLL = " << hPos.NLL << ", vertex = " << vPos[0]<<", "<<vPos[1]<<", "<<vPos[2]<<", "<<vPos[3]<<std::endl;
#endif
#endif
	}
	
	tBestReconstructedVertexPosition = tReconstructedVertexPosition[0]; 
			
#ifdef CHECK_TO_TRUE_VTX
	double distanceToTrue = GetDistance(tBestReconstructedVertexPosition,fTrueVtxPos);
	
#ifdef VERBOSE_MIN
	std::cout << "Fine candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << tBestReconstructedVertexPosition[VTX_T] - fTrueVtxPos[VTX_T] << "ns, with a min NLL = " <<minNLL << std::endl;
	std::cout << "Time of the true event = " << fTrueVtxPos[VTX_T] << std::endl;
#endif
#endif
	clock_t timeEnd=clock();
	
	//if(verbose) 
		std::cout<<"MinimizeVertex: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<std::endl<<std::endl;

	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
}


std::vector< std::vector<double> > BQFitter::MinimizeVertex_Main(	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
		
		
	std::vector<std::thread> lThreadList;
	std::vector< std::vector<double> > lOutputFinal;
	
	int iCand_Step = nCandidates / N_THREAD;
	
	if ( iCand_Step < 1 ) iCand_Step = 1;
	
	for ( int iStart = 0; iStart < nCandidates; iStart += iCand_Step ) {
			
		std::thread tThrd( MinimizeVertex_CallThread, iStart,iCand_Step, 
							initialVertex, limits, stepSize, nhits, nCandidates, 
							tolerance, verbose, likelihood, average, lowerLimit, upperLimit, directionality);
							
		lThreadList.push_back(std::move(tThrd));
	}
	
	
	for ( unsigned int iIdx=0; iIdx < lThreadList.size(); iIdx++ ) {
		// Wait thread
		if ( lThreadList[iIdx].joinable() )
			lThreadList[iIdx].join();
	}
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	lOutputFinal = fThreadOutput;
	fThreadOutput.clear();
	mtx.unlock();
		
	std::sort(lOutputFinal.begin(), lOutputFinal.end(), SortOutputVector); 	

	while((int) lOutputFinal.size() > tolerance) {	
		lOutputFinal.pop_back();
	}
	
	
	return lOutputFinal;
}
		
//here is the function where minimization gonna take place. The loop is inside this function
void BQFitter::MinimizeVertex_thread(	
		int iStart, int iIte,	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool /*likelihood*/,bool /*average*/,double lowerLimit, double upperLimit, int directionality){
		
	int iEnd = (iIte+iStart);
	if ( iEnd > nCandidates ) iEnd = nCandidates;

	std::vector<struct FitPosition>  tVtxContainer;

	TFitter * minimizer = new TFitter(8);//4=nb de params?
	TMinuit * minuit = minimizer->GetMinuit();
	
	double arglist[20];
	int err=0;
	//arglist[0]=0;
	double p1=verbose-1;
	minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
	minuit->SetErrorDef(1);
	minuit->SetPrintLevel(-1); // quiet mode
	arglist[0]=2;
	minuit->mnexcm("SET STR",arglist,1,err);//set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
	minimizer->SetFCN(MinuitLikelihood);//here is function to minimize

	double Mig[2]={1e6,1e0};//maxcalls and tolerance
	
	minimizer->SetParameter(0,"vertex0",0,stepSize		  ,-limits[0]			,limits[0]		);
	minimizer->SetParameter(1,"vertex1",0,stepSize		  ,-limits[1]			,limits[1]		);
	minimizer->SetParameter(2,"vertex2",0,stepSize		  ,-limits[2]			,limits[2]		);
	minimizer->SetParameter(3,"vertex3",0,stepSize/fLightSpeed,-limits[3]/fLightSpeed	,limits[3]/fLightSpeed	);
	
	minimizer->SetParameter(4,"PMTConfiguration",0,0,0,0); // Useless now
	minimizer->SetParameter(5,"nhits",nhits,nhits,nhits-1,nhits+1);
	minimizer->SetParameter(6,"lowerLimit",lowerLimit,5e-1,-7,-2);
	minimizer->SetParameter(7,"upperLimit",upperLimit,5e-1,2,10);
	minimizer->SetParameter(8,"expoSigma",100,5,0,1e3);
	minimizer->SetParameter(9,"directionality",directionality,1,0,2);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);
	

	for(int icand=iStart;icand<iEnd;icand++){
						
		minimizer->SetParameter(0,"vertex0",initialVertex[icand][0],stepSize		  ,initialVertex[icand][0]-limits[0]			,initialVertex[icand][0]+limits[0]		);
		minimizer->SetParameter(1,"vertex1",initialVertex[icand][1],stepSize		  ,initialVertex[icand][1]-limits[1]			,initialVertex[icand][1]+limits[1]		);
		minimizer->SetParameter(2,"vertex2",initialVertex[icand][2],stepSize		  ,initialVertex[icand][2]-limits[2]			,initialVertex[icand][2]+limits[2]		);
		minimizer->SetParameter(3,"vertex3",initialVertex[icand][3],stepSize/fLightSpeed  ,initialVertex[icand][3]-limits[3]/fLightSpeed	,initialVertex[icand][3]+limits[3]/fLightSpeed	);
	
		minimizer->SetParameter(4,"PMTConfiguration",0,0,0,0); // Useless now
		minimizer->SetParameter(5,"nhits",nhits,nhits,nhits-1,nhits+1);
		minimizer->SetParameter(6,"lowerLimit",lowerLimit,5e-1,-7,-2);
		minimizer->SetParameter(7,"upperLimit",upperLimit,5e-1,2,10);
		minimizer->SetParameter(8,"expoSigma",100,5,0,1e3);
		minimizer->SetParameter(9,"directionality",directionality,1,0,2);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);
		
		minimizer->ExecuteCommand("MIGRAD",Mig,2);

		std::vector<double> tCandRecoVtxPos(4,0.);
		
		tCandRecoVtxPos[0] = minimizer->GetParameter(0);
		tCandRecoVtxPos[1] = minimizer->GetParameter(1);
		tCandRecoVtxPos[2] = minimizer->GetParameter(2);
		tCandRecoVtxPos[3] = minimizer->GetParameter(3);
		
		struct FitPosition hPos;
		hPos.Vtx = tCandRecoVtxPos;
		hPos.NLL = 0;
		
		// Get minuit
		minuit=minimizer->GetMinuit();
		double fmin, fedm, errdef;
		int npar,nparx,istat;
		minuit->mnstat(fmin,fedm,errdef,npar,nparx,istat);
		
		hPos.NLL = fmin;
		
		tVtxContainer.push_back(hPos);
	}	 
	 
	std::sort(tVtxContainer.begin(),tVtxContainer.end(),SortingNLL());
	while((int) tVtxContainer.size() > tolerance) {	
		tVtxContainer.pop_back();
	}
	
	// Check if mutex is already lock, if not -> lock
	// We don't want simultaneous modification of the common output
	mtx.lock();
	for ( unsigned int iPos = 0; iPos < tVtxContainer.size(); iPos++ ) {
		std::vector<double> vPos(5,0.);
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


void BQFitter::ComputeInTimeHit(int iType) {
	// Get n hit in a given time window (after tof correction)
	
	int iNbr = fHitInfo.size();
	std::vector<PMTHitExt> tHitTmp;
	
	bool bCond = (iType != AllPMT);
	
	// Compute time of flight and some geometric information
	for ( int iHit=0; iHit < iNbr; iHit++ ) {
		int    iPMT    = fHitInfo[iHit].PMT;
		
		
		// Skip Hybrid PMT is the option is set
		if ( bCond ) {
			if ( GetPMTType(iPMT) != iType ) continue; 
		}
		PMTHitExt tNewHit;
		
		tNewHit.PMT  = iPMT;
		tNewHit.T    = fHitInfo[iHit].T;
		tNewHit.Q    = fHitInfo[iHit].Q;
		
		tNewHit.Dist  = GetDistance( fPMT_Info[iPMT], fRecoVtxPosFinal[0] );
		tNewHit.ToF   = tNewHit.T - (tNewHit.Dist / CNS2CM);
		
		tNewHit.NormX = (fPMT_Info[iPMT][0] - fRecoVtxPosFinal[0][0])/tNewHit.Dist;
		tNewHit.NormY = (fPMT_Info[iPMT][1] - fRecoVtxPosFinal[0][1])/tNewHit.Dist;
		tNewHit.NormZ = (fPMT_Info[iPMT][2] - fRecoVtxPosFinal[0][2])/tNewHit.Dist;
		
		// Not needed
		tNewHit.theta = 0.;
		tNewHit.phi   = 0.;
		
		// Store hit
		tHitTmp.push_back(tNewHit);
	}
	
	iNbr = tHitTmp.size();
	
	//std::cout << " InTime computation " << iNbr << " type " << iType << std::endl;
	
	// Clear in time Hit list 
	fInTime20.clear();
	fInTime30.clear();
	fInTime50.clear();

	// Make tmp in time hit list
	std::vector< PMTHitExt > tInTime20;
	std::vector< PMTHitExt > tInTime30;
	std::vector< PMTHitExt > tInTime50;
	
	std::sort(tHitTmp.begin(),tHitTmp.end(),SortingToF());
	
	for ( int iHit=0; iHit < iNbr; iHit++ ) {
	
		PMTHitExt hHit = tHitTmp[iHit];
			
		double dTime = hHit.ToF;
		
		tInTime20 .push_back(hHit);
		tInTime30 .push_back(hHit);
		tInTime50 .push_back(hHit);
		
		while ( (dTime - tInTime20 [0].ToF) > 20. ) {
			tInTime20 .erase(tInTime20 .begin());
		}
		while ( (dTime - tInTime30 [0].ToF) > 30. ) {
			tInTime30 .erase(tInTime30 .begin());
		}
		while ( (dTime - tInTime50 [0].ToF) > 50. ) {
			tInTime50 .erase(tInTime50 .begin());
		}
		
		if ( tInTime20 .size() > fInTime20 .size() )
			fInTime20  = tInTime20;
			
		if ( tInTime30 .size() > fInTime30 .size() )
			fInTime30  = tInTime30;
			
		if ( tInTime50 .size() > fInTime50 .size() )
			fInTime50  = tInTime50;
	}
	
	//std::cout << " InTime computation " << fInTime20.size() << " " << fInTime30.size() << " " << fInTime50.size() << std::endl;
}

void BQFitter::ComputeDistanceFromWall() {
	
	double dR = fTankRadius - fRecoVtx_R;
	double dZ = fTankHalfHeight - abs(fRecoVtx_Z);
	
	fRecoVtx_Wall = std::min(dR,dZ);
}

double BQFitter::GetPMTAngleEfficiency(double dCosPM) {
	// From sklib/coseffsk.F

	if ( dCosPM > 1. ) {
		std::cout << "WARNING: Abnormal dCosPM " << dCosPM << std::endl;
		dCosPM = 1.;	
	}
	if ( dCosPM < 0. ) {
		dCosPM = 0.;	
	}
	
	//double sk1_a[4] = { 0.3542135, 0.5104711, 0.1160429, 0.0118482 }; // SK-I
	double sk2_a[4] = { 0.205349 , 0.523981 , 0.389951 ,-0.131959  }; // SK-II, SK-III, SK-IV, SK-V
	
	double * data = sk2_a;
	
	return (  data[0] 
		+ data[1] * dCosPM
		+ data[2] * dCosPM * dCosPM
		+ data[2] * dCosPM * dCosPM * dCosPM );
}

std::vector<double> BQFitter::TransformInVector(std::vector<double> lDir, double dPhi, double dCosTheta ) {

	std::vector<double> lOut(3,0.);
	
	double dSinTheta = sqrt( 1. - (dCosTheta*dCosTheta) );
	
	double dX = dSinTheta * cos(dPhi);
	double dY = dSinTheta * sin(dPhi);
	double dZ = dCosTheta;
	
	if ( !(abs(lDir[0]) < 1e-2 && abs(lDir[1]) < 1e-2) ) {
	
		double dRDir = sqrt(  	lDir[0] * lDir[0] + 
					lDir[1] * lDir[1] + 
					lDir[2] * lDir[2] );
					
		double dCosA = lDir[2] / dRDir;
		double dSinA = sqrt( 1. - (dCosTheta*dCosTheta) );
		
		double dCosB = lDir[0] / dRDir / dSinA;
		double dSinB = lDir[1] / dRDir / dSinA;
		
		double dXp =  dX * dCosA + dZ * dSinA;
		double dYp =  dY;
		double dZp = -dX * dSinA + dZ * dCosA;
		
		lOut[0] = dXp * dCosB - dYp * dSinB;
		lOut[1] = dXp * dSinB + dYp * dCosB;
		lOut[2] = dZp;
		
		return lOut;
	}
	
	if ( lDir[2] >= 0 ) {
		lOut[0] = dX;
		lOut[1] = dY;
		lOut[2] = dZ;
	}
	else {
		lOut[0] = -dX;
		lOut[1] = -dY;
		lOut[2] = -dZ;
	}
	
	return lOut;
}


void BQFitter::LoadDirConstant() {

	// From lowe/sklowe/dirlkflf.F
	fEP10[ 0] =  959.1;
	fEP10[ 1] = 1077.29;
	fEP10[ 2] = 1186.07;
	fEP10[ 3] = 1305.84;
	fEP10[ 4] = 1437.70;
	fEP10[ 5] = 1582.88;
	fEP10[ 6] = 1742.72;
	fEP10[ 7] = 2000.  ;
	fEP10[ 8] = 1846.76;
	fEP10[ 9] = 1709.92;
	fEP10[10] = 1584.08;
	fEP10[11] = 1469.25;
	fEP10[12] = 1365.42;
	fEP10[13] = 1272.60;
	fEP10[14] = 1190.78;
	fEP10[15] = 1119.97;
	fEP10[16] = 1060.17;
	fEP10[17] = 1011.37;
	fEP10[18] =  973.57;
	fEP10[19] =  946.78;
	fEP10[20] =  931.  ;
		
	// From lowe/sklowe/lfneweff_sk3_final.F
	fEffHitLambda[ 0] =  7000.;
	fEffHitLambda[ 1] =  8000.;
	fEffHitLambda[ 2] =  9000.;
	fEffHitLambda[ 3] = 10000.;
	fEffHitLambda[ 4] = 11000.;
	fEffHitLambda[ 5] = 12000.;
	fEffHitLambda[ 6] = 13000.;
	fEffHitLambda[ 7] = 14000.;
	fEffHitLambda[ 8] = 15000.;
	fEffHitLambda[ 9] = 16000.;
	fEffHitLambda[10] = 17000.;
	fEffHitLambda[11] = 18000.;
	
	
	// From lowe/sklowe/lfoccor_3.F
	// Occupancy correction table
	// Likely need to be modified to include mPMT later
	int iFullOccupancy  	= 3;
	fNearPMT		= 10;
	for ( int i=0; i < fNearPMT; i++ ) {
		for ( int j=0; j < fNearPMT; j++ ) {
		
			double d1 = (double) i;
			double d2 = (double) j;
			if ( i == j )		
				fOccupancyTable[i][j] = (double) iFullOccupancy;
			else if ( j > i || i == 0 || j == 0 ) 
				fOccupancyTable[i][j] = 0;
			else
				fOccupancyTable[i][j] = - d1 * log(1.0 - d2 / d1) / d2;
				// value: - N/m * log(1 - m/N)
				// 		N: number of total alive PMT
				//		m: number of hit PMT
		}
	}
	
	
}

double BQFitter::GetDirKS(std::vector<double> lDir) {

	double dCosTh = lDir[2];
	//double dSinTh = sqrt( (1. - dCosTh) * (1. * dCosTh) );
	double dNorm  = sqrt( lDir[0] * lDir[0] + lDir[1] * lDir[1] );
	double dCosPh = lDir[0] / dNorm;
	double dSinPh = lDir[1] / dNorm;
	
	int n50 = fInTime50.size();
	
	std::vector<double> dPhi;
	
	for ( int iHit=0; iHit < n50; iHit++ ) {
	
		double dXX = dSinPh * fInTime50[iHit].NormX - dCosPh * fInTime50[iHit].NormY;
		double dYY = 	  dCosTh * dCosPh * fInTime50[iHit].NormX
				+ dCosTh * dSinPh * fInTime50[iHit].NormY
				-          dSinPh * fInTime50[iHit].NormZ;
				
		dPhi.push_back( atan2( dYY, dXX ) * 180. / 3.1416 + 180. );
	}
	
	std::sort(dPhi.begin(),dPhi.end());
	
	double dNormDeg =  360. / (double) n50;
	double dKSMax   = -360.;
	double dKSMin   =  360.;
	
	for ( int iHit=0; iHit < n50; iHit++ ) {
	
		double dKS = dPhi[iHit] - dNormDeg * (iHit+1);
		
		if ( dKS >= dKSMax ) dKSMax = dKS;
		if ( dKS <= dKSMin ) dKSMin = dKS;
	}
	
	double dDirKS = (dKSMax-dKSMin) / 360.;
	
	return dDirKS;
}

double BQFitter::GetDir_lkflf(double dCos, int iE, bool bSimple) {

	double dOutput = 0.;
	
	if ( bSimple ) {
		// From lowe/sklowe/dirlkflf.F
		
		if ( dCos < -1 || dCos > 1. ) {
			std::cout << "WARNING: dCos out of range in GetDir_lkflf simple " << dCos << std::endl;
			dCos = dCos < 0. ? -1. : 1.; 
		}
		
		if ( dCos <= 0.2 ) {
			dOutput = exp( 5.61 + 1.36 * dCos ) / 2000.;
		}
		else if ( dCos <= 0.6 ) {
			dOutput = exp( 5.42 + 2.41 * dCos ) / 2000.;
		}
		else {
			int    iCos  = (int) ( dCos * 50. ) - 30;
			double dFrac = ( dCos * 50. - 30. ) - (double) iCos;
			
			
			if ( iCos == 20 ) 
				dOutput = fEP10[20];
			else 
				dOutput = fEP10[iCos] + dFrac * ( fEP10[iCos+1] - fEP10[iCos] );
				
			dOutput /= 2000.;
		}
		
		return dOutput;
	}
	
	// From low/sklowe/dirlkflf_ene_sk4.F
	int iCos = (int) ( ( dCos + 1.0 ) * 50 );
	return fHitPat[iE][iCos-1];
}

std::vector<double> BQFitter::GetDirection(double dEnergy, std::vector<PMTHitExt> tInTime, bool bSimple, double dCutOff) {
	// Compute direction and direction prob
	//  Input needed (for dir4): energy, n20, false, 0.2
	//  Input needed (for dir2): energy, n30, true,  0.01
	 
	std::vector<double> tOutput(4,0);
	int nHit = tInTime.size();
	
	int iE = 0;
	
	// Get vector index for GetDir_lkflf
	if ( !bSimple ) {
		if      ( dEnergy <  3.5 ) iE =  0;
		else if ( dEnergy <  4.5 ) iE =  1;
		else if ( dEnergy <  5.5 ) iE =  2;
		else if ( dEnergy <  6.5 ) iE =  3;
		else if ( dEnergy <  7.5 ) iE =  4;
		else if ( dEnergy <  8.5 ) iE =  5;
		else if ( dEnergy <  9.5 ) iE =  6;
		else if ( dEnergy < 10.5 ) iE =  7;
		else if ( dEnergy < 11.5 ) iE =  8;
		else if ( dEnergy < 12.5 ) iE =  9;
		else if ( dEnergy < 13.5 ) iE = 10;
		else if ( dEnergy < 14.5 ) iE = 11;
		else if ( dEnergy < 15.5 ) iE = 12;
		else if ( dEnergy < 16.5 ) iE = 13;
		else if ( dEnergy < 19.0 ) iE = 14;
		else if ( dEnergy < 21.0 ) iE = 15;
		else if ( dEnergy < 23.0 ) iE = 16;
		else if ( dEnergy < 25.0 ) iE = 17;
		else if ( dEnergy < 27.0 ) iE = 18;
		else if ( dEnergy < 40.0 ) iE = 19;
		else if ( dEnergy < 60.0 ) iE = 20;
		else			   iE = 21;
	}
	
	
	if ( nHit > 1 ) {
	
		// Variables:		
		double dDeltaXsum = 0;
		double dDeltaYsum = 0;
		double dDeltaZsum = 0;
	
		for ( int iHit=0; iHit < nHit; iHit++ ) {			
			dDeltaXsum += tInTime[iHit].NormX;
			dDeltaYsum += tInTime[iHit].NormY;
			dDeltaZsum += tInTime[iHit].NormZ;
		}
		
		double dDist = sqrt( 	  dDeltaXsum * dDeltaXsum 
					+ dDeltaYsum * dDeltaYsum 
					+ dDeltaZsum * dDeltaZsum );
					
		tOutput[0] = dDeltaXsum / dDist;
		tOutput[1] = dDeltaYsum / dDist;
		tOutput[2] = dDeltaZsum / dDist;
	
		double dProb = 0;
		
		for ( int iHit=0; iHit < nHit; iHit++ ) {
			int iPMT = tInTime[iHit].PMT;
			
			double dCosTh =   tOutput[0] * tInTime[iHit].NormX 
					+ tOutput[1] * tInTime[iHit].NormY 
					+ tOutput[2] * tInTime[iHit].NormZ;
			// 20080718 tentative fix by y.takeuchi  (is this ok??)
			if ( dCosTh > 1. && dCosTh < 1.0001 ) {
				dCosTh = 1.;
			}
			// end
			
			double dTmp = GetDir_lkflf(dCosTh,iE,bSimple);
			if ( dTmp <= dCutOff ) dTmp = dCutOff;
			
			
			double dCosPM = - fPMT_Info[iPMT][3] * tInTime[iHit].NormX
					- fPMT_Info[iPMT][4] * tInTime[iHit].NormY
					- fPMT_Info[iPMT][5] * tInTime[iHit].NormZ;
					
			//double dRR    = GetDistance(fPMT_Info[iPMT],fRecoVtxPosFinal[0]);
			
			// 20080718 tentative fix by y.takeuchi  (is this ok??)
			if ( dCosPM > 1. && dCosPM < 1.0001 ) {
				dCosPM = 1.;
			}
			// end
			
			double dCorEff = dCosPM / GetPMTAngleEfficiency(dCosPM) * tInTime[iHit].Q;
			
			dProb += (log10(dTmp) * dCorEff);
		}
		// Save prob
		tOutput[3] = dProb;
		
		double dTheta[4] = 	{
						20. *3.1416/180.,
						 9. *3.1416/180.,
						 4. *3.1416/180.,
						 1.6*3.1416/180.,
					};
					
		for ( int iStep = 0; iStep < 4; iStep++ ) {
			
			double dCosTheta = cos(dTheta[iStep]);
			
			int iTry  = 0;
			int iMove = 1;
			
			while ( iTry < 9 && iMove == 1 ) {
				iMove  = 0;
				iTry  += 1;
				
				for ( int iPhi=0; iPhi < 6; iPhi++ ) {
					
					double dPhi = float(iPhi) / 6. * 2. * 3.1416;
					
					std::vector<double> lWDir = TransformInVector(tOutput,dPhi,dCosTheta);
					
					double dProb2 = 0.;
					
					for ( int iHit=0; iHit < nHit; iHit++ ) {
						int iPMT = tInTime[iHit].PMT;
						
						double dCosTh =   lWDir[0] * tInTime[iHit].NormX
								+ lWDir[1] * tInTime[iHit].NormY 
								+ lWDir[2] * tInTime[iHit].NormZ;
								
						// 20080718 tentative fix by y.takeuchi  (is this ok??)
						if ( dCosTh > 1. && dCosTh < 1.0001 ) {
							dCosTh = 1.;
						}
						// end
						
						double dTmp = GetDir_lkflf(dCosTh,iE,bSimple);
						if ( dTmp <= dCutOff ) dTmp = dCutOff;
												
						double dCosPM = - fPMT_Info[iPMT][3] * tInTime[iHit].NormX
								- fPMT_Info[iPMT][4] * tInTime[iHit].NormY
								- fPMT_Info[iPMT][5] * tInTime[iHit].NormZ;
								
						//double dRR    = GetDistance(fPMT_Info[iPMT],fRecoVtxPosFinal[0]);
						
						// 20080718 tentative fix by y.takeuchi  (is this ok??)
						if ( dCosPM > 1. && dCosPM < 1.0001 ) {
							dCosPM = 1.;
						}
						// end			
						
						double dCorEff = dCosPM / GetPMTAngleEfficiency(dCosPM) * tInTime[iHit].Q;
						
						dProb2 += (log10(dTmp) * dCorEff);
					}
					
					if ( dProb2 > dProb ) {
						dProb   = dProb2;
						iMove   = 1;
						tOutput = lWDir;
						tOutput.push_back( dProb );
					}	
				} // for iPhi
			} // while
		} // for iStep
	}
	
	return tOutput;
}

void BQFitter::MakeAnalysis(int iType) {

	// Compute In time hit (n50, n30, ...)
	this->ComputeInTimeHit(iType);
			
	fVtxReco_dirKS = 2.;
			
	fVtxReco_dirKS 		= 0.;
	fRecoVtx_E_Simple	= 0.;
	fRecoVtx_E		= 0.;
	fVtxReco_dir.clear();
	fVtxReco_dirSimple.clear();
	
			
	if ( fRecoVtx_Wall > 0. && fInTime50.size() > 0 ) {
			
		// Compute simple direction
		// original: 
		//	call lfdir2(bsvertex, bsdir_lfdir2, bsdirks)
		// qmin = .2, qmax = 25., cutoff = .01, twindow = 30.
		fVtxReco_dirSimple = this->GetDirection(0., fInTime30, true, 0.01);
	
		// Compute dirKS from Simple direction
		if ( fVtxReco_dirSimple[3] != -1 ) {
			// Check goodness of direction
			// original:
			//	call lfdirks(bsvertex, bsdir, bsdirks)
			fVtxReco_dirKS = this->GetDirKS(fVtxReco_dirSimple);			
		}				
	}
	else {
		fVtxReco_dirSimple.assign(4,0.);
	}
}

struct BQFitter::FitterOutput BQFitter::MakeFit(bool bHybrid) {

	// Get Hits total
	int iHitsTotal = fHitInfo.size();
	
	int iPMTConfiguration = 0; // Normal PMT
	if ( bHybrid ) {
		iPMTConfiguration = 1; // Hybrid configuration
	}
	
	//Proceed to the fit.
	double * tLimits = new double[4];
	
	struct FitterOutput fOutput;
	
	fOutput.Vtx[0]      = 0.;
	fOutput.Vtx[1]      = 0.;
	fOutput.Vtx[2] 	    = 0.;
	fOutput.Time        = 0.;
	fOutput.Good        = 0.;
	
	fOutput.InTime      = 0;
	
	fOutput.True_NLLDiff  = 0;
	fOutput.True_TimeDiff = 0;
	fOutput.True_TistDiff = 0;
	
	fLastLowerLimit       = 0.;
	fLastUpperLimit       = 0.;
	
	
	fRecoVtx_X = 0; 
	fRecoVtx_Y = 0;  
	fRecoVtx_Z = 0;  
	fRecoVtx_T = 0;  
	fRecoVtx_R2 = 0; 
	fRecoVtx_R = 0;
	
	fRecoVtx_Wall = 0;
	
	fInTime20.clear();
	fInTime30.clear();
	fInTime50.clear();
			
	/*
	std::cout << " nMult PMT = " << fWCGeo->GetWCNumPMT(true) << std::endl;
	for(int a=0;a<iHitsTotal;a++) {
		std::cout<<"Hit [ " << a << " ] PMT = " <<fHitInfo[a].PMT  << " T = " << fHitInfo[a].T  << " Q = " << fHitInfo[a].Q << std::endl;
		std::cout<< " PMT: (" 
					<< fPMT_Info[fHitInfo[a].PMT][0] << " , " 
					<< fPMT_Info[fHitInfo[a].PMT][1] << " , " 
					<< fPMT_Info[fHitInfo[a].PMT][2] << " )( " 
					<< fPMT_Info[fHitInfo[a].PMT][3] << " , " 
					<< fPMT_Info[fHitInfo[a].PMT][4] << " , " 
					<< fPMT_Info[fHitInfo[a].PMT][5] << " ) " 
					<< std::endl;
	}
	  */    	
	if ( true ) {
	
		if(iHitsTotal>0){
		
			TStopwatch timer;
			timer.Reset();
			timer.Start();
				
			//2.a. Coarse vertex search from scratch
			double dStepSize=fSearchVtxStep;//in centimeters
			int iTolerance=60;//30
			
			//std::vector< std::vector<double> > tRecoVtxPos = this->SearchVertex(iHitsTotal,iTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false/*fUseDirectionality*/);
			std::vector< std::vector<double> > tRecoVtxPos = this->SearchVertex_Main(iHitsTotal,iTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false/*fUseDirectionality*/);
			
						
			if(VERBOSE==2){
				for(unsigned int a=0;a<tRecoVtxPos.size();a++) std::cout<<"Candidate vertex [ " << a << " ] = " <<tRecoVtxPos[a][0] << " , " << tRecoVtxPos[a][1] << " , " << tRecoVtxPos[a][2] << " , " << tRecoVtxPos[a][3] << std::endl;
			}
			timer.Stop();
			//std::cout << "SearchVertex took: " << timer.RealTime() << std::endl;
			////////////////////////////////////////
#ifdef KILLFAKE
			bool bTrueIncluded=false;
			for(int a=0;a<iTolerance;a++){
	    
				double dTime2True=0;
				double dDist2True=0;
				
				for(int j=0;j<4;j++){
					if(j==3) dTime2True = tRecoVtxPos[a][j]-fTrueVtxPos[j];
					else{
						//std::cout<<"Distance in "<<j<<" is ="<<tRecoVtxPos[a][j]<<","<<fTrueVtxPos[j]<<std::endl;
						//std::cout<<"Distance in "<<j<<" is ="<<pow(tRecoVtxPos[a][j]-fTrueVtxPos[j],2)<<std::endl;
						dDist2True+=(tRecoVtxPos[a][j]-fTrueVtxPos[j])*(tRecoVtxPos[a][j]-fTrueVtxPos[j]);
					}
				}
				dDist2True=TMath::Sqrt(dDist2True);
				if(dDist2True < 2*dStepSize && dTime2True < 2*dStepSize/fLightSpeed) bTrueIncluded=true;
			}
			if(!bTrueIncluded) continue;
#endif
		
			double dStepSizeFinal=10;
			int iToleranceFinal=1;
			if(fStepByStep){
				//2.b. Refined vertex search around candidate found previously
				
				//double dStepSizeFine = 300;
				//int iToleranceFine=40;
				//for(int i=0;i<4;i++) tLimits[i+1] = dStepSize;//particleStart[i]*1e-2;//onversion to meter
		    
				//double ** tRecoVtxPosFine = this->SearchVertexFine(tRecoVtxPos,tLimits,dStepSizeFine,iHitsTotal,iTolerance,iToleranceFine,VERBOSE);
				//double ** tRecoVtxPosFine = this->SearchVertexFine(tRecoVtxPos,tLimits,dStepSizeFine,iHitsTotal,iTolerance,iToleranceFine,VERBOSE,true,false,-100,500);
				//if(VERBOSE == 2){
					//for(int a=0;a<toleranceFine;a++) std::cout<<"Fine candidate vertex time = "<<tRecoVtxPosFine[a][0]<<std::endl;
				//}
				////////////////////////////////////////
		    
				//2.c. More refined vertex search around candidate found previously
				double dStepSizeFine2=100;
				int iToleranceFine2=20;
				for(int i=0;i<4;i++) tLimits[i] = dStepSize;//particleStart[i]*1e-2;//onversion to meter
				std::vector< std::vector<double> > tRecoVtxPosFine2 = this->SearchVertexFine(tRecoVtxPos,tLimits,dStepSizeFine2,iHitsTotal,iTolerance,iToleranceFine2,VERBOSE,true,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,fUseDirectionality);
				//double ** tRecoVtxPosFine2 = this->SearchVertexFine(tRecoVtxPosFine,tLimits,dStepSizeFine2,iHitsTotal,iToleranceFine,iToleranceFine2,1,true,false,-100,500);
				//////////////////////////////////////////
		    
				//2.d. Final refined vertex search around candidate found previously
				double dStepSizeFine3=20;
				int iToleranceFine3=1;
				for(int i=0;i<4;i++) tLimits[i] = dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
				std::vector< std::vector<double> > tRecoVtxPosFine3 = this->SearchVertexFine(tRecoVtxPosFine2,tLimits,dStepSizeFine3,iHitsTotal,iToleranceFine2,iToleranceFine3,VERBOSE,true,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,fUseDirectionality);
				//double ** tRecoVtxPosFine3 = this->SearchVertexFine(tRecoVtxPosFine2,tLimits,dStepSizeFine3,iHitsTotal,iToleranceFine2,iToleranceFine3,VERBOSE,true,false,-100,500);
				//for(int i=0;i<4;i++) tLimits[i] = 0;//dStepSizeFine;//particleStart[i]*1e-2;//onversion to meter
				//double ** tRecoVtxPosFineTest = this->SearchVertexFine(tRecoVtxPosFine3,tLimits,dStepSizeFine3,iHitsTotal,1,toleranceFine3,2,true,false);
				if(VERBOSE == 1){
					for(int a=0; a<iToleranceFine3; a++) std::cout<<"Final candidate vertex time = "<<tRecoVtxPosFine3[a][0]<<", x = "<<tRecoVtxPosFine3[a][1]<<", y = "<<tRecoVtxPosFine3[a][2]<<", z="<<tRecoVtxPosFine3[a][3]<<std::endl;
				}
				//////////////////////////////////////////	
				fPDFNorm_fullTimeWindow = this->SplineIntegral(fSplineTimePDFQueue[iPMTConfiguration],fSTimePDFLimitsQueueNegative_fullTimeWindow,fSTimePDFLimitsQueuePositive_fullTimeWindow);
				for(int i=0;i<4;i++) tLimits[i] = 100;//dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
				//double ** fRecoVtxPosFinal = tRecoVtxPosFine3;//this->SearchVertexFine(,tLimits,dStepSizeFinal,iHitsTotal,iToleranceFine3,iToleranceFinal,VERBOSE,true,false);
				//double ** fRecoVtxPosFinal = minimizeVertex(fTrueVtxPosDouble,tLimits,dStepSizeFinal,iHitsTotal,iToleranceFine3,iToleranceFinal,VERBOSE,true,false,-100,100);
				fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPosFine3,tLimits,dStepSizeFinal,iHitsTotal,iToleranceFine3,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,fUseDirectionality);
				if(VERBOSE == 1){
					for(int a=0;a<iToleranceFinal;a++){
						std::cout<<"Final candidate vertex time = "<<fRecoVtxPosFinal[a][0]<<", x = "<<fRecoVtxPosFinal[a][1]<<", y = "<<fRecoVtxPosFinal[a][2]<<", z="<<fRecoVtxPosFinal[a][3]<<std::endl;
					}
				}
			}
			else {
			
				timer.Reset();
				timer.Start();
				fPDFNorm_fullTimeWindow = this->SplineIntegral(fSplineTimePDFQueue[iPMTConfiguration],fSTimePDFLimitsQueueNegative_fullTimeWindow,fSTimePDFLimitsQueuePositive_fullTimeWindow);
				for(int i=0;i<4;i++) tLimits[i] = 2*dStepSize;
				
				//fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPos,tLimits,dStepSizeFinal,iHitsTotal,iTolerance,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,fUseDirectionality);
				fRecoVtxPosFinal = this->MinimizeVertex_Main(tRecoVtxPos,tLimits,dStepSizeFinal,iHitsTotal,iTolerance,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,fUseDirectionality);
				
				//for(int a=0;a<iToleranceFinal;a++) std::cout<<"Final candidate vertex time = "<<fRecoVtxPosFinal[a][3]<<", x = "<<fRecoVtxPosFinal[a][0]<<", y = "<<fRecoVtxPosFinal[a][1]<<", z="<<fRecoVtxPosFinal[a][2]<<std::endl;
				//double dStepSizeFine = 50;
				//int iToleranceFine=30;
				//double ** tRecoVtxPosFine = this->MinimizeVertex(tRecoVtxPos,tLimits,dStepSizeFine,iHitsTotal,iTolerance,iToleranceFine,VERBOSE,true,false,-100,500);
				//for(int i=0;i<4;i++) tLimits[i] = dStepSizeFine*4;
				//fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPosFine,tLimits,dStepSizeFinal,iHitsTotal,iTolerance,iToleranceFinal,VERBOSE,true,false,-100,500);
				//if(VERBOSE == 1){
				//	for(int a=0;a<iToleranceFinal;a++) std::cout<<"Final candidate vertex time = "<<fRecoVtxPosFinal[a][0]<<", x = "<<fRecoVtxPosFinal[a][1]<<", y = "<<fRecoVtxPosFinal[a][2]<<", z="<<fRecoVtxPosFinal[a][3]<<std::endl;
				//}
				timer.Stop();
				//std::cout << "Minimizer took: " << timer.RealTime() << std::endl;
			}

		/*
			//Bonus:
			//double dStepSizeFinal=20;
			//int iToleranceFine3=1;
			for(int i=0;i<4;i++) tLimits[i] = 0;//dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
			double ** fRecoVtxPosFinal;
			//for(int a=0;a<iToleranceFine3;a++){
				fRecoVtxPosFinal = this->SearchVertexFine(tRecoVtxPosFine3,tLimits,dStepSizeFine3,iHitsTotal,iToleranceFine3,iToleranceFine3,1,true,-1,2);
				if(VERBOSE == 1){
					for(int a=0;a<iToleranceFine3;a++){
						//std::cout<<"Updated final candidate vertex time = "<<fRecoVtxPosFinal[a][0]<<", x = "<<fRecoVtxPosFinal[a][1]<<", y = "<<fRecoVtxPosFinal[a][2]<<", z="<<fRecoVtxPosFinal[a][3]<<std::endl;
					}
				}
			//}
			//
		*/

			//3.a. Test the NLL of the true vertex to compare with reco
			
			//if(VERBOSE) std::cout<<"Test of the NLL at the candidate vertex point"<<std::endl;
			//for(int i=0;i<4;i++) tLimits[i] = 0;//dStepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
			//double ** reconstructedVertexPositionTest0 = this->SearchVertexFine(&(tRecoVtxPosFine[0]),tLimits,dStepSizeFineFine,iHitsTotal,1,100,VERBOSE,true);
			
#ifdef CHECK_TO_TRUE_VTX		
			if(VERBOSE) std::cout<<"Test of the NLL at the true vertex point"<<std::endl;
			for(int i=0;i<4;i++) tLimits[i] = 0;//dStepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
			std::vector< std::vector<double> > tTrueVtxPosFinal = this->SearchVertexFine(fTrueVtxPosDouble,tLimits,dStepSize,iHitsTotal,1,1,VERBOSE,true);
			if(VERBOSE) std::cout<<"done"<<std::endl;
#ifdef OUTPUT_TREE					
			//DistanceToTrueXNLL[iPMTConfiguration]->Fill(dDist2True);
#endif
			/////////////////////////////////////////////////////////////

		  
			// Bonus: Probe the NLL space near the true vertex:
			if(bFirstEvent){
				for(int i=0;i<4;i++) tLimits[i] = 2000;//dStepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
				double dStepSizeTrue=500;
				int iToleranceTrue=200;
				std::cout<<"Scan"<<std::endl;
				//std::vector< std::vector<double> > reconstructedVertexPosition = searchVertex(fTankRadius,fTankHeight,fIntegrationTimeWindow,dStepSize,iHitsTotal,iTolerance,VERBOSE,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false/*fUseDirectionality*/);
				//std::vector< std::vector<double> > tTrueVtxPosScan = searchVertex(fTrueVtxPosDouble,tLimits,dStepSizeTrue,iHitsTotal,1,iToleranceTrue,2,true,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,2/*fUseDirectionality*/);
				std::vector< std::vector<double> > tTrueVtxPosScan = this->SearchVertexFine(fTrueVtxPosDouble,tLimits,dStepSizeTrue,iHitsTotal,1,iToleranceTrue,2,true,false,-10,20,2/*fUseDirectionality*/);


				for(int a=0;a<iToleranceTrue;a++){
					double dTime2True=0;
					double nllDifference=0;
					double dDist2True=0;
					for(int j=0;j<4;j++){
						if(j==3) {
							dTime2True = tTrueVtxPosScan[a][j]-fTrueVtxPos[j];
						}
						else{
							//std::cout<<"Distance in "<<j<<" is ="<<tRecoVtxPosFine3[a][j]<<","<<fTrueVtxPos[j]<<std::endl;
							//std::cout<<"Distance in "<<j<<" is ="<<pow(tRecoVtxPosFine3[a][j]-fTrueVtxPos[j],2)<<std::endl;
							dDist2True+=(tTrueVtxPosScan[a][j]-fTrueVtxPos[j])*(tTrueVtxPosScan[a][j]-fTrueVtxPos[j]);
						}
					}
					dDist2True=TMath::Sqrt(dDist2True);
					
					nllDifference=tTrueVtxPosScan[a][4]-tTrueVtxPosFinal[0][4];
					
					fOutput.True_NLLDiff  = nllDifference;
					fOutput.True_TimeDiff = dTime2True;
					fOutput.True_TistDiff = dDist2True;
#ifdef OUTPUT_TREE					
					fitNLLprofile_currentEvent[iPMTConfiguration]->Fill(dTime2True,dDist2True,nllDifference);
					normfitNLLprofile_currentEvent[iPMTConfiguration]->Fill(dTime2True,dDist2True,1);
#endif
					//dDist2True=TMath::Sqrt(dDist2True)*1e2;//Conversion back to cm
				}
			}
			///////////////////////////////////////////////
		
		  
			//3.b. For final candidate, test how far it is from true vertex
			double dTime2True=0;
			double dDist2True=0;
			double nllDifference=0;

			if(VERBOSE) std::cout<<"done"<<std::endl;
			for(int a=0;a<iToleranceFinal;a++){
				if(VERBOSE) std::cout<<"Final candidate vertex time = "<<fRecoVtxPosFinal[a][0]<<std::endl;
				if(a==0){
					for(int j=0;j<4;j++){
						if(j==3) dTime2True = fRecoVtxPosFinal[a][j]-fTrueVtxPos[j];
						else{
							if(VERBOSE) std::cout<<"Distance in "<<j<<" is ="<<fRecoVtxPosFinal[a][j]<<","<<fTrueVtxPos[j]<<std::endl;
							//std::cout<<"Distance in "<<j<<" is ="<<pow(tRecoVtxPosFine3[a][j]-fTrueVtxPos[j],2)<<std::endl;
							dDist2True+=(fRecoVtxPosFinal[a][j]-fTrueVtxPos[j])*(fRecoVtxPosFinal[a][j]-fTrueVtxPos[j]);
						}
					}
					nllDifference=fRecoVtxPosFinal[a][4]-tTrueVtxPosFinal[0][4];
					//dDist2True=TMath::Sqrt(dDist2True)*1e2;//Conversion back to cm
				}
			}
			
			if(VERBOSE) std::cout<<"done"<<std::endl;
			dDist2True=TMath::Sqrt(dDist2True);
			if(tTrueVtxPosFinal[0][4]!=0) nllDifference/=tTrueVtxPosFinal[0][4];
			if(VERBOSE) std::cout<<"Vertex found has a time difference of "<<dTime2True<<"ns and a distance of "<<dDist2True<<"cm from the true vertex, and nll difference ="<<nllDifference<<std::endl;
			if(dDist2True>500) {
				std::cout<<"Vertex found has a time difference of "<<dTime2True<<"ns and a distance of "<<dDist2True<<"cm from the true vertex, and nll difference ="<<nllDifference<<std::endl;
				std::cout<<"code 666"<<std::endl;
			}
#ifdef OUTPUT_TREE					
			fitTimeToTrue[iPMTConfiguration]->Fill(dTime2True);
			fitDistanceToTrue[iPMTConfiguration]->Fill(dDist2True);
			fit4DdistanceXNLLdifference[iPMTConfiguration]->Fill(dDist2True,nllDifference);
			std::cout << std::endl;
#endif
			
			//////////////////////////////////////////////////////////////////
#ifdef OUTPUT_TREE					
			for(int j=0;j<4;j++){
				/* if(j<3){
					bsvertex[j] = fRecoVtxPosFinal[0][j+1];
					mcvertex[j] = tTrueVtxPosFinal[0][j+1];
				}
				else{
					bsvertex[3] = fRecoVtxPosFinal[0][0];
					mcvertex[3] = tTrueVtxPosFinal[0][0];
				}*/
				if(j==3){
					bsvertex[0] = fRecoVtxPosFinal[0][3];
					mcvertex[0] = fTrueVtxPos[3];
				}
				else if(j==0){
					bsvertex[3] = fRecoVtxPosFinal[0][0];
					mcvertex[3] = fTrueVtxPos[0];
				}
				else{
					bsvertex[j] = fRecoVtxPosFinal[0][j];
					mcvertex[j] = fTrueVtxPos[j];
				}
		    
			}
#endif

#endif

/*
			std::cout<<"Vertex found has a time difference of "<< fTrueVtxPos[3]-fRecoVtxPosFinal[0][3] <<"ns and a distance of "<<GetDistance(fRecoVtxPosFinal[0],fTrueVtxPos)<<"cm from the true vertex, and nll ="<<fRecoVtxPosFinal[0][4]<<std::endl;
			std::cout<<"Reco Vtx (x,y,z,t) " << fRecoVtxPosFinal[0][0] << ", " << fRecoVtxPosFinal[0][1] << ", " << fRecoVtxPosFinal[0][2] << ", " << fRecoVtxPosFinal[0][3] << std::endl;
			std::cout<<"True Vtx (x,y,z,t) " << fTrueVtxPos[0] << ", " << fTrueVtxPos[1] << ", " << fTrueVtxPos[2] << ", " << fTrueVtxPos[3] << std::endl;
*/			

						
			int iInTime = 0;
			
			for(int ihit = 0; ihit < iHitsTotal; ihit++){
				int iPMT = fHitInfo[ihit].PMT;
				double hitTime = fHitInfo[ihit].T;
				double distance = GetDistance(fPMT_Info[iPMT],fRecoVtxPosFinal[0]);
				double tof = distance / fLightSpeed;
				double residual = hitTime - tof - fRecoVtxPosFinal[0][3];
				if(residual > fSTimePDFLimitsQueueNegative && residual < fSTimePDFLimitsQueuePositive){
					iInTime++;
				}
			}
			
			// Fill variables for additionnal analysis
			fRecoVtx_X	= fRecoVtxPosFinal[0][0]; 
			fRecoVtx_Y	= fRecoVtxPosFinal[0][1];  
			fRecoVtx_Z	= fRecoVtxPosFinal[0][2];  
			fRecoVtx_T	= fRecoVtxPosFinal[0][3];  
			fRecoVtx_R2	= fRecoVtx_X * fRecoVtx_X + fRecoVtx_Y * fRecoVtx_Y; 
			fRecoVtx_R	= sqrt(fRecoVtx_R2);
			fRecoVtx_Good	= fRecoVtxPosFinal[0][4];  
			
			// Fill Output
			fOutput.Vtx[0] 	= fRecoVtx_X;
			fOutput.Vtx[1] 	= fRecoVtx_Y;
			fOutput.Vtx[2] 	= fRecoVtx_Z;
			fOutput.Time   	= fRecoVtx_T;
			fOutput.InTime 	= iInTime;
			fOutput.Good 	= fRecoVtx_Good;
			
			
			//std::cout << " Final vertex " << fOutput.Vtx[0] << " " << fOutput.Vtx[1] << " " << fOutput.Vtx[2] << " Time " << fOutput.Time << std::endl;
			
			this->ComputeDistanceFromWall();
			fOutput.Wall  = fRecoVtx_Wall;
			
			this->MakeAnalysis(NormalPMT);
			
			fOutput.n50 	    [NormalPMT]    = fInTime50.size();
			fOutput.dirKS 	    [NormalPMT]    = fVtxReco_dirKS;
			fOutput.dir 	    [NormalPMT][0] = fVtxReco_dirSimple[0];
			fOutput.dir 	    [NormalPMT][1] = fVtxReco_dirSimple[1];
			fOutput.dir 	    [NormalPMT][2] = fVtxReco_dirSimple[2];
			fOutput.dir_goodness[NormalPMT]    = fVtxReco_dirSimple[3];
			
			this->MakeAnalysis(MiniPMT);
			
			fOutput.n50 	    [MiniPMT  ]    = fInTime50.size();
			fOutput.dirKS 	    [MiniPMT  ]    = fVtxReco_dirKS;
			fOutput.dir 	    [MiniPMT  ][0] = fVtxReco_dirSimple[0];
			fOutput.dir 	    [MiniPMT  ][1] = fVtxReco_dirSimple[1];
			fOutput.dir 	    [MiniPMT  ][2] = fVtxReco_dirSimple[2];
			fOutput.dir_goodness[MiniPMT  ]    = fVtxReco_dirSimple[3];
			
			this->MakeAnalysis(AllPMT);
			
			fOutput.n50 	    [AllPMT   ]    = fInTime50.size();
			fOutput.dirKS 	    [AllPMT   ]    = fVtxReco_dirKS;
			fOutput.dir 	    [AllPMT   ][0] = fVtxReco_dirSimple[0];
			fOutput.dir 	    [AllPMT   ][1] = fVtxReco_dirSimple[1];
			fOutput.dir 	    [AllPMT   ][2] = fVtxReco_dirSimple[2];
			fOutput.dir_goodness[AllPMT   ]    = fVtxReco_dirSimple[3];
			
#ifdef OUTPUT_TREE								
			bstree->Fill();
#endif
		}
	}
	
	return fOutput;
}
