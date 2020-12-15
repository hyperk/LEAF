/*****************************************************************************************************/
/**	LEAF.cc											**/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)					**/
/**	Original author: Benjamin Quilain								**/
/**	Date: December 18th 2019									**/
/**	Desc: Low-E Fitter for Hyper-K								**/
/*****************************************************************************************************/

#include "LEAF.hh"
#include "TStopwatch.h"

double LEAF::fSTimePDFLimitsQueueNegative = 0;
double LEAF::fSTimePDFLimitsQueuePositive = 0;
LEAF* LEAF::myFitter=NULL;

std::mutex mtx;

/*****************************************************************************************************/

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
	
	NLL=LEAF::GetME()->FindNLL_Likelihood(vertexPosition,nhits,lowerLimit,upperLimit,true,false,directionality);	
	
	//timer.Stop();7
	//std::cout << " Likelihood: " << timer.RealTime() << std::endl;

}


void MinimizeVertex_CallThread(	
		int iStart, int iIte,	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
	
	LEAF::GetME()->MinimizeVertex_thread(iStart, iIte,	
		initialVertex, limits, stepSize, nhits,
		nCandidates, tolerance, verbose, likelihood, average, lowerLimit,  upperLimit, directionality);
	
}	

void SearchVertex_CallThread(
			int iStart, int iIte,
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
			
	LEAF::GetME()->SearchVertex_thread(iStart, iIte,	
		nhits, tolerance, likelihood, lowerLimit, upperLimit, directionality);
			
}
/*****************************************************************************************************/

LEAF::LEAF() {

	myFitter = this;
	fThread  = N_THREAD;
}

LEAF::~LEAF() {
	
	//fTrueVtxPosDouble.clear();
	delete fRand;
	
}

LEAF* LEAF::GetME() {

	if ( myFitter ) return myFitter;

	myFitter = new LEAF();
	return myFitter;
}

void LEAF::DeleteME() {

	if ( myFitter ) delete myFitter;
	
}

/*****************************************************************************************************/
void LEAF::Initialize(const Geometry* lGeometry) {
	
	fGeometry = lGeometry;
	
	// Get ID PMT info list
	fPMTList = fGeometry->GetPMTList();
	
	// Set DarkNoise:
	fDarkRate_ns[NormalPMT] 	= fGeometry->pmt_dark_rate[HKAA::kIDPMT_BnL];
	fDarkRate_ns[MiniPMT]		= fGeometry->pmt_dark_rate[HKAA::kIDPMT_3inch];
	
	// Set Tank size:
	fTankRadius 			= fGeometry->detector_radius;
	fTankHeight			= fGeometry->detector_length;
	fTankHalfHeight		= fTankHeight / 2.;
	
	this->Init();
}
/*****************************************************************************************************/

void LEAF::Init() {

	// Parameters
	fStepByStep 					= false; // Step by step mode was not test and is not supported by multithreading
	fUseDirectionality 				= false;
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
	
	fMinimizeLimitsNegative			= -700;
	fMinimizeLimitsPositive			= 1000;
	
	fSearchVtxStep					= 300;//in cm, the step size for coarse grid search
	fSearchVtxTolerance 				= 60;//Number of candidate vertex that are kept after coarse grid search	
	fHitTimeLimitsNegative				= -5;//-10;//-5;//-5
	fHitTimeLimitsPositive				= 7;//15;//7;//7
	
	fIntegrationTimeWindow 			= 50;//in ns
	
	// Compute light speed:
	float fCVacuum 				= 3e8*1e2 / 1e9;//speed of light, in centimeter per ns.
	float fNIndex 					= 1.373;//1.385;//1.373;//refraction index of water
	fLightSpeed 					= fCVacuum/fNIndex;
	
	/*
	fTrueVtxPos.clear();
	fTrueVtxPosDouble.clear();
	fTrueVtxPos.resize(5,0.);
	fTrueVtxPosDouble.push_back(fTrueVtxPos);
	*/
	
	fPDFNorm_fullTimeWindow 			= 0.;
	
	//fPositionSize				= 0;
			
	// Get random generator
	fRand  					= new TRandom3();
	
	this->LoadSplines();

	// Make position lists
      	this->MakePositionList();
}



void LEAF::LoadSplines() {

	TFile *fSplines, *fSplines2;
	if(fHighEnergy){
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_HE.root","read");//To generate with my code ProduceWSPlots.c
		fSplines2 = fSplines;
	}
	else{
		fSplines = new TFile("${LEAFDIR}/inputs/timePDF_DRnew_Large.root","read");//To generate with my code ProduceWSPlots.c
		fSplines2 = new TFile("${LEAFDIR}/inputs/timePDF_Directionality_DRnew.root","read");//To generate with my code ProduceWSPlots.c
	}
	
	// Prevent TGraph2D to be append to TFile (this is needed as we are doing multiple copy of TGraph2D)
	// For a strange reason sometimes the TGraph2D is see as a TH1 and stored in an open TFile
	TH1::AddDirectory(false);

	for(int pmtType=0;pmtType<NPMT_CONFIGURATION;pmtType++){
		//Load 1D t-tof splines
		fSplineTimePDFQueue[pmtType]    = (TSpline3*) fSplines->Get(Form("splineExpoQueue%d_%d",0,pmtType));
		fSplineTimePDFDarkRate[pmtType] = (TSpline3*) fSplines->Get(Form("splineDR%d_%d",0,pmtType));

		fSTimePDFLimitsQueueNegative_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmin();
		fSTimePDFLimitsQueuePositive_fullTimeWindow = fSplineTimePDFQueue[pmtType]->GetXmax();
		//std::cout<<"Min="<<fSTimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<fSTimePDFLimitsQueuePositive_fullTimeWindow<<std::endl;
		if(fSTimePDFLimitsQueueNegative < fSTimePDFLimitsQueueNegative_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmin();//To avoid to use PDF where it is not defined
		if(fSTimePDFLimitsQueuePositive > fSTimePDFLimitsQueuePositive_fullTimeWindow) fSTimePDFLimitsQueuePositive = fSplineTimePDFQueue[pmtType]->GetXmax();//same here.

		//Load 3D directionality histograms.
		for(int pmtGroup=0;pmtGroup<HKAA::kmPMT_Groups;pmtGroup++){
			gPMTDirectionality_2D[pmtType][pmtGroup] = (TGraph2D*) fSplines2->Get(Form("gPMTDirectionality_2D_%d_%d_%d",0,pmtType,pmtGroup));
		}
		fDistResponsePMT[pmtType] = (TF1*) fSplines2->Get(Form("fDistResponsePMT_pmtType%d",pmtType));
	}
}

/*
void LEAF::SetTrueVertexInfo(std::vector<double> vtx, double time) {

	fTrueVtxPosDouble.clear();
	fTrueVtxPos.clear();
	
	fTrueVtxPos.resize(5,0.);
	
	fTrueVtxPos[0] = vtx[0];
	fTrueVtxPos[1] = vtx[1];
	fTrueVtxPos[2] = vtx[2];
	fTrueVtxPos[3] = time;
	fTrueVtxPos[4] = 0.;
	
	fTrueVtxPosDouble.push_back(fTrueVtxPos);
}*/
/*****************************************************************************************************/
// Spline Integral functions:
double LEAF::SplineIntegral(TSpline3 * s,double start,double end,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*s->Eval(i+stepSize/2);
	}
	return integral;
}
double LEAF::SplineIntegralAndSubstract(TSpline3 * s0,TSpline3 * s1,double start,double end,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*(s0->Eval(i+stepSize/2)-s1->Eval(i+stepSize/2));
	}
	return integral;
}

double LEAF::SplineIntegralExpo(TSpline3 * s,double start,double end,double sigma,double stepSize){
	double integral=0;
	for(double i=start;i<end;i+=stepSize){
		integral+=stepSize*s->Eval(i+stepSize/2)*TMath::Gaus(i+stepSize/2,0,sigma);
	}
	return integral;
}

/*****************************************************************************************************/

void LEAF::VectorVertexPMT(std::vector<double> vertex, int iPMT, double* dAngles) {
		
	// Guillaume 2020/05/20:
	// mPMT referencial is computed once for all PMT in LoadPMTInfo() and MakeMPMTReferencial(iPMT)
	// Guillaume 2020/11/20:
	// Reference is moved to outside LEAF (HKManager or hk-AstroAnalysis)
  	
	PMTInfo lPMTInfo = (*fPMTList)[iPMT];
  	int iPMTTop = lPMTInfo.mPMT_RefTube;
	PMTInfo lPMTInfoTop = (*fPMTList)[iPMTTop];
  
	//5. Now we have our referential, we should just calculate the angles of the PMT to vertex position vector in this referential.
	//a. calculate the PMT to vertex position vector.
	
	double dVtx_PMTRef[3];
	
	dVtx_PMTRef[0] = vertex[0] - lPMTInfoTop.Position[0];
	dVtx_PMTRef[1] = vertex[1] - lPMTInfoTop.Position[1];
	dVtx_PMTRef[2] = vertex[2] - lPMTInfoTop.Position[2];
	
	double dLengthVtx = Astro_GetLength(dVtx_PMTRef);
	GeoTools::Normalize(dVtx_PMTRef);
		
	if( VERBOSE >= 3 ){
		std::cout << "Vertex position = " << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << std::endl;
		std::cout << "Vertex position from PMT = " << dVtx_PMTRef[0] << ", " << dVtx_PMTRef[1] << ", " << dVtx_PMTRef[2] << std::endl;
	}
	
	//b. Then extract Theta and Phi:
	double dCosTheta= Astro_GetScalarProd(dVtx_PMTRef,lPMTInfo.mPMT_RefZ);
	double dTheta	= TMath::ACos(dCosTheta);
	
	double dPhi 	= 0.;
	
	if( (*fPMTList)[iPMT].mPMT_TubeNum == HKAA::kmPMT_TopID ) {
		//Phi is not defined in that case..
	}
	else {
		//We know x=cosPhi x sinTheta and y=sinPhi x sinTheta
		double dX 	= Astro_GetScalarProd(dVtx_PMTRef,lPMTInfo.mPMT_RefX);
		double dY 	= Astro_GetScalarProd(dVtx_PMTRef,lPMTInfo.mPMT_RefY);
		double dTanPhi	= dY/dX;
		dPhi		= TMath::ATan(dTanPhi);
		
		//tan is symetric from -pi/2 to +pi/2
		if( dX == 0 ){
			if( dY < 0 )	dPhi = -TMath::Pi()/2.;
			else 		dPhi =  TMath::Pi()/2.;
		}
		if( dX < 0 ) 		dPhi += TMath::Pi(); //With this, angle become defines between -pi/2 to -pi/2 +2pi.
		
		//We wish to bring this from 0 to 2pi:
		if( dPhi < 0 ) 		dPhi += 2*TMath::Pi();
		
		//Actually, we have a symmetry in Phi between [0,pi] and [pi,2pi]. So, we will just define Phi in [0,pi] modulo pi
		if( dPhi > TMath::Pi() )dPhi = TMath::Pi() - (dPhi - TMath::Pi());
	}
	
	dAngles[0] = dPhi   * 180./TMath::Pi();
	dAngles[1] = dTheta * 180./TMath::Pi();
	dAngles[2] = dLengthVtx;
	dAngles[3] = 1; //calculateWeight(dLengthVtx,PMTradius[pmtType],Theta,verbose);
    
	if( VERBOSE >= 3 ){
		std::cout << "Angles Phi = " << dAngles[0] << ", Theta = " << dAngles[1] << std::endl;
	}
	
}
double LEAF::FindDirectionTheta(std::vector<double> vertex,int tubeNumber, int verbose) {

	//clock_t timeStart=clock();
	
	PMTInfo lPMTInfo = (*fPMTList)[tubeNumber];
	
	int pmt_number_in_mpmt = lPMTInfo.mPMT_TubeNum;
	double particleRelativePMTpos[3];
	for(int j=0;j<3;j++) particleRelativePMTpos[j] = lPMTInfo.Position[j] - vertex[j];
	if(verbose == 3){
		std::cout<<"Original PMT number = "<<pmt_number_in_mpmt<<std::endl;  
		std::cout<<"Direction = "<<lPMTInfo.Orientation[0]<<", "<<lPMTInfo.Orientation[1]<<", "<<lPMTInfo.Orientation[2]<<std::endl;
	}

	double norm=0;
	double ps=0;
	for(int j =0;j<3;j++){
		ps   += (-particleRelativePMTpos[j]) * lPMTInfo.Orientation[j];
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

void LEAF::MakeEventInfo(double lowerLimit, double upperLimit) {

	if ( fLastLowerLimit == lowerLimit && fLastUpperLimit == upperLimit ) return;

	fLastLowerLimit = lowerLimit;
	fLastUpperLimit = upperLimit;
	
	for(int i=0; i < NPMT_CONFIGURATION; i++) {
		fEventInfo[i].hits 		= 0;
		fEventInfo[i].SignaloverNoise	= 0.;
		fEventInfo[i].NoiseIntegral	= 0.;
		fEventInfo[i].SignalIntegral 	= 0.;
	}
	
	int iHitTotal = fHitCollection->Size();
	
	for(int iHit=0; iHit < iHitTotal; iHit++ ) {
		//Hit lHit = fHitInfo[iHit];		
		Hit lHit = fHitCollection->At(iHit);
		
		int iPMT = lHit.PMT;
		int iType = Astro_GetPMTType(iPMT);
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
		
		if(VERBOSE>=3) std::cout<<"nhits="<<fEventInfo[i].hits<<", DR average="<<signalDR
					<<", signal over noise="<<fEventInfo[i].SignaloverNoise
					<<", signal integral="  <<fEventInfo[i].SignalIntegral
					<<", DR integral="	<<fEventInfo[i].NoiseIntegral
					<<",in integral="	<<fEventInfo[i].SignalIntegral/fEventInfo[i].NoiseIntegral<<std::endl;
	
		
		//if(VERBOSE>=2){
		//	std::cout<<"Lower limit = "<<lowerLimit<<", proba = "<<fSplineTimePDFQueue[i]->Eval(lowerLimit)<<std::endl;
		//	std::cout<<"Upper limit = "<<upperLimit<<", proba = "<<fSplineTimePDFQueue[i]->Eval(upperLimit)<<std::endl;
		//}
	
	}
	

}

void LEAF::MakePositionList() {
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

/*****************************************************************************************************/

double LEAF::FindNLL_Likelihood(std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality){
	
	//TStopwatch timer;
	//timer.Reset();
	//timer.Start();
	
	double NLL = 0;
		
	for(int ihit = 0; ihit < nhits; ihit++){
		//Hit lHit = fHitInfo[ihit];		
		Hit lHit = fHitCollection->At(ihit);
		
		int iPMT = lHit.PMT;
		double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		
		int pmtType = Astro_GetPMTType(iPMT);
		double distance = Astro_GetDistance(lPMTInfo.Position,vertexPosition); 
		
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
		
		condition = residual > lowerLimit && residual < upperLimit;
		
		if(condition){
		
			proba = fSplineTimePDFQueue[pmtType]->Eval(residual);
			
#ifdef VERBOSE_NLL
			if(VERBOSE>=3){
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
				if(VERBOSE>=3) std::cout<<"proba after scaling="<<proba<<std::endl;
#endif
			}
			
		}
		else if(killEdges) {
			if ( residual >= upperLimit ){
				proba=fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Upper limit = "<<upperLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else if( residual <= lowerLimit ){
				//continue;
				proba=fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Lower limit = "<<lowerLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
		}
		else{
			proba=0;
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
		
		NLLdir=this->FindNLLDirectionality(vertexPosition, nhits, VERBOSE,lowerLimit,upperLimit);
		
#ifdef VERBOSE_NLL
		if(VERBOSE>=2) std::cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdir)<<", dir 2 = "<<NLLdir<<std::endl;
#endif
		if(directionality == 1) NLL+=NLLdir;
		else if(directionality == 2) NLL=NLLdir;
	}
	
	//timer.Stop();
	
	//std::cout << " FindNLL_Likelihood full: " << timer.RealTime() << std::endl;

	return NLL;
}


double LEAF::FindNLL_NoLikelihood(std::vector<double> vertexPosition, int nhits, double /*lowerLimit*/, double /*upperLimit*/, bool /*killEdges*/, bool /*scaleDR*/, int /*directionality*/){
	double NLL = 0;
				
	//std::cout << " Find NLL " << nhits << std::endl;
	//std::cout << " First Hit " << fHitCollection->At(0).PMT << " " << vertexPosition.size() << std::endl;	
	for(int ihit = 0; ihit < nhits; ihit++){
		//Hit lHit = fHitInfo[ihit];		
		Hit lHit = fHitCollection->At(ihit);
	
		//std::cout << " NLL Hit " << ihit << " " << lHit.PMT << std::endl;
		int iPMT = lHit.PMT;		
	
		double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		
		int pmtType = Astro_GetPMTType(iPMT);
		double distance = Astro_GetDistance(lPMTInfo.Position,vertexPosition);
		
		double tof = distance / fLightSpeed;
		double residual = hitTime - tof - vertexPosition[3];
#ifdef KILLHALF
		if(pmtType==1){
			double t=fRand->Uniform(0,1);
			if(t>0.5) continue;
		}
#endif

		bool bCondition = (residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive) || (pmtType == 1 && !fLimit_mPMT);

      		//std::cout << " HIT " << ihit << " Has " <<   bCondition << " " << pmtType << " " << fLimit_mPMT << " " << residual << " ( " << hitTime << " - " << tof << " - " << vertexPosition[3] << " ) "  << distance << " / " << fLightSpeed << " " << lPMTInfo.Position[0] << " "<< lPMTInfo.Position[1] << " "<< lPMTInfo.Position[2] << " " << iPMT << std::endl;
      
		if( bCondition ){
			NLL++;
#ifdef VERBOSE_NLL
			if(VERBOSE>=3){
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
	
	
	//std::cout<<"NLL="<<NLL<<std::endl;

	return NLL;
}


double LEAF::FindNLL(std::vector<double> vertexPosition, int nhits, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality){
	double NLL = 0;
	
	TStopwatch timer;
	timer.Reset();
	timer.Start();
	
	for(int ihit = 0; ihit < nhits; ihit++){
		//Hit lHit = fHitInfo[ihit];		
		Hit lHit = fHitCollection->At(ihit);
		
		int iPMT = lHit.PMT;
		double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		
		int pmtType = Astro_GetPMTType(iPMT);
		double distance = Astro_GetDistance(lPMTInfo.Position,vertexPosition);
		
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
			
			condition = residual > lowerLimit && residual < upperLimit;
			if(condition){
				proba = fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3){
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
					if(VERBOSE>=3) std::cout<<"proba after scaling="<<proba<<std::endl;
				}
			}
			else if(killEdges && residual >= upperLimit){
				proba=fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Upper limit = "<<upperLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else if(killEdges && residual <= lowerLimit){
				proba=fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3 && pmtType==1) std::cout<<"PMT type = "<< pmtType <<", Lower limit = "<<lowerLimit<<", residual = "<<residual<<", proba = "<<proba<<std::endl;
#endif
			}
			else{
				proba=0;
			}
		}
		else{
			bool condition;
			if(pmtType == 0 || (pmtType == 1 && fLimit_mPMT) ) condition = residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive;
			else condition = true;
			if(condition){
				NLL++;//= fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if(VERBOSE>=3){
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
	//if(VERBOSE>=2) cout<<"NLL="<<NLL<<endl;
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
		if(VERBOSE>=2) std::cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdir)<<", dir 2 = "<<NLLdir2<<std::endl;
		if(directionality == 1) NLL+=NLLdir2;
		else if(directionality == 2) NLL=NLLdir2;
	}
	timer.Stop();
	std::cout<<"FindNLL: Total time = "<< timer.RealTime() << std::endl;
	timer.Reset();
	

	return NLL;
}

double LEAF::FindNLLDirectionality(std::vector<double> vVtxPos, int nhits, int verbose, double /*lowerLimit*/, double /*upperLimit*/) {

	double NLL = 0;
	std::vector<double> vDirection(2,0.);//Return phi and theta.
	
	// Interpolate isn't compatible with multi-thread, need to use mutex which lead to long deadtime.
	// Copy TGraph2D 
	mtx.lock();
	static thread_local TGraph2D gPMTDirectionality_2D_local_0 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][0]);
	static thread_local TGraph2D gPMTDirectionality_2D_local_1 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][1]);
	static thread_local TGraph2D gPMTDirectionality_2D_local_2 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][2]);
	mtx.unlock();
	
	TGraph2D* tDirectionality[3];
	tDirectionality[0] = &gPMTDirectionality_2D_local_0;
	tDirectionality[1] = &gPMTDirectionality_2D_local_1;
	tDirectionality[2] = &gPMTDirectionality_2D_local_2;
	
		
	for(int ihit = 0; ihit < nhits; ihit++){
		//Hit lHit = fHitInfo[ihit];		
		Hit lHit = fHitCollection->At(ihit);
		
		int iPMT = lHit.PMT;
		
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		int pmtType = Astro_GetPMTType(iPMT);

		if(pmtType==0) continue;
		bool condition = true;
		
		if(condition){
			double vPMTVtx[4];
			this->VectorVertexPMT(vVtxPos,iPMT,vPMTVtx);
			
			double dPhi   = vPMTVtx[0];
			double dTheta = vPMTVtx[1];
			double dDist  = vPMTVtx[2];
			
			int pmtGroup = lPMTInfo.mPMT_Group;
			
			double proba = 0;
			proba = tDirectionality[pmtGroup]->Interpolate(dPhi,dTheta);
						
			double dDistCorr = fDistResponsePMT[pmtType]->Eval(dDist);
			
			proba *= dDistCorr;
      
			if(proba==0){
				proba=1e-20;
				if(verbose) std::cout<<"We are at proba = 0, theta = "<< dTheta <<", proba used = "<< proba <<std::endl;
			}
			NLL += -TMath::Log(proba);
      
			if(VERBOSE>=3){
				//std::cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << vDirection[1] << ", proba="<<proba<<std::endl;
				std::cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << dTheta << ", proba="<<proba<<std::endl;
			}
		}
	}
	if(VERBOSE>=2) std::cout<<"NLL directionnel="<<NLL<<std::endl;
	return NLL;
}
/*****************************************************************************************************/

//Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
//The step between each grid points is given by stepSize (in centimeters).
//We give the list of hit information through the 2D array fHitInfo to make the code faster, instead of using the whole code information
//The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower). 
//The number of output vertices can be chosen using "tolerance".
std::vector< std::vector<double> > LEAF::SearchVertex(int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	
	std::vector<double> tBestReconstructedVertexPosition;
	
	//double ** reconstructedVertexPosition = new double[4];
		
	std::vector<struct FitPosition>  tVtxContainer;
	
	clock_t timeStart=clock();
		
	if ( fPositionList.size() == 0 ) std::cout << " Position list not filled " << std::endl;
	if(VERBOSE>=3) std::cout<<"Number of hits for coarse grid search = "<<nhits<<std::endl;
	for ( unsigned int iPos = 0; iPos < fPositionList.size(); iPos++ ) {
		
  		struct FitPosition hPos;
  		hPos.Vtx = fPositionList[iPos];
  		hPos.NLL = 0;
  							
		if ( likelihood ) 	hPos.NLL = this->FindNLL_Likelihood  (hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		else 			hPos.NLL = this->FindNLL_NoLikelihood(hPos.Vtx,nhits,lowerLimit,upperLimit,false,false,directionality);
		
		hPos.NLL += fRand->Uniform(0,1e-5);//Here only not to have exactly the same value for 2NLL, and be removed from map
	  
	  	tVtxContainer.push_back(hPos);
	  	
		if(VERBOSE>=3) std::cout << "Test Pos: " << hPos.Vtx[3] << " " << hPos.Vtx[0] << " " << hPos.Vtx[1] << " " << hPos.Vtx[2] << " NLL = " << hPos.NLL << " " << likelihood << std::endl;
					
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
		std::cout<<"NLL = " << hPos.NLL << ", vertex = " << vPos[0]<<", "<<vPos[1]<<", "<<vPos[2]<<", "<<vPos[3]<<std::endl;
#endif
	}
	
	clock_t timeEnd=clock();
	
	//if(VERBOSE) 
		std::cout<<"SearchVertex: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<< " " << fPositionList.size() << std::endl;
	
	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
}

std::vector< std::vector<double> > LEAF::SearchVertex_Main(	
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
		
		
	std::vector<std::thread> lThreadList;
	std::vector< std::vector<double> > lOutputFinal;
	
	int iCand_Step = (int) fPositionList.size() / fThread;
	
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
void LEAF::SearchVertex_thread(
			int iStart, int iIte,
			int nhits,int tolerance,bool likelihood,double lowerLimit, double upperLimit,int directionality){
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	
	std::vector<struct FitPosition>  tVtxContainer;
			
	if ( fPositionList.size() == 0 ) std::cout << " Position list not filled " << std::endl;
	if(VERBOSE>=3) std::cout<<"Number of hits for coarse grid search = "<<nhits<<std::endl;	
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
		//hPos.NLL += fRand->Uniform(0,1e-5);//Here only not to have exactly the same value for 2NLL, and be removed from map
	  
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
std::vector< std::vector<double> > LEAF::SearchVertexFine(std::vector< std::vector<double> > initialVertex,double * limits, double stepSize,int nhits,int nCandidates,int tolerance,int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
	
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
			if(VERBOSE>=2) std::cout << "Time = " << time << "ns" << std::endl;
			for(double x = initialVertex[icand][VTX_X] - limits[1]; x <= initialVertex[icand][VTX_X] + limits[1];x+=stepSize){
				if(VERBOSE>=2) std::cout << "x = " << x << "m" << std::endl;
				for(double y = initialVertex[icand][VTX_Y] - limits[2]; y <= initialVertex[icand][VTX_Y] + limits[2];y+=stepSize){
					if(VERBOSE>=2) std::cout << "y = " << y << "m" << std::endl;
					for(double z = initialVertex[icand][VTX_Z] - limits[3]; z <= initialVertex[icand][VTX_Z] + limits[3];z+=stepSize){
						if(VERBOSE>=2) std::cout << "z = " << z << "m" << std::endl;

						std::vector<double> tVtxPosition(4,0.);//In centimeters
						tVtxPosition[0] = x;//radius*TMath::Cos(angle);
						tVtxPosition[1] = y;//radius*TMath::Sin(angle);
						tVtxPosition[2] = z;//height;
						tVtxPosition[3] = time;

						double NLL=0;
						if(average){
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
								NLL+=nll;
								//vNLL[irand]=nll;
							}
	      
							if(fAveraging!=0) NLL/=fAveraging;
						}
						else NLL=this->FindNLL(tVtxPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);

	    
						std::vector<double> candidateReconstructedVertexPosition(4,0.);
						for(int i=0;i<4;i++) candidateReconstructedVertexPosition[i] = tVtxPosition[i];
	    
						vertexContainer[NLL] = candidateReconstructedVertexPosition;
						if( (int) vertexContainer.size() > tolerance){
							vertexContainer.erase(--vertexContainer.rbegin().base());
						}
						
						if(NLL<minNLL){
							minNLL=NLL;
							//std::cout << "minNLL = "<<minNLL << std::endl;
							tBestReconstructedVertexPosition = tVtxPosition;
						}
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
			}
		counter++;
	}
	
	clock_t timeEnd=clock();
	if(verbose) 
		std::cout<<"SearchVertexFine: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<std::endl<<std::endl;
  
	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
  
}


//here is the function where minimization gonna take place. The loop is inside this function
std::vector< std::vector<double> > LEAF::MinimizeVertex(std::vector< std::vector<double> > initialVertex,double * limits, double stepSize,int nhits,int nCandidates,int tolerance,int verbose,bool /*likelihood*/,bool /*average*/,double lowerLimit, double upperLimit, int directionality){
	if(VERBOSE>=2) std::cout<<"Minimizer"<<std::endl;
	//2. How to set the search?
	//double stepSize = 0.5;//in m
	//double * reconstructedVertexPosition = new double[4];
	std::vector<double> tBestReconstructedVertexPosition;
	clock_t timeStart=clock();

	std::vector<struct FitPosition>  tVtxContainer;

	if(VERBOSE>=2) std::cout<<"Minimizer"<<std::endl;
	TFitter * minimizer = new TFitter(8);//4=nb de params?
	TMinuit * minuit = minimizer->GetMinuit();
	
	double arglist[20];
	int err=0;
	//arglist[0]=0;
	double p1=verbose-1;
	minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
	minuit->SetErrorDef(1);
	//minuit->SetPrintLevel(-1); // quiet mode
	if(VERBOSE<2) {
		minuit->mnexcm("SET NOWarnings",0,0,err);
	}
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
	minimizer->SetParameter(10,"thread",0,0,0,0);
	minimizer->FixParameter(4);
	minimizer->FixParameter(5);
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);
	minimizer->FixParameter(10);

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
		//minimizer->SetParameter(8,"scaling factor",1,1,0,1e3);
		//minimizer->SetParameter(8,"signal pe",1,1,0,1e3);
		minimizer->SetParameter(9,"directionality",directionality,1,0,2);
		minimizer->SetParameter(10,"thread",0,0,0,0);
		minimizer->FixParameter(4);
		minimizer->FixParameter(5);
		minimizer->FixParameter(6);
		minimizer->FixParameter(7);
		minimizer->FixParameter(8);
		minimizer->FixParameter(9);
		minimizer->FixParameter(10);
      
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
		
#ifdef VERBOSE_MIN
		if(verbose>=1) std::cout<<"NLL="<<hPos.NLL<<std::endl;
#endif
		tVtxContainer.push_back(hPos);

		if(VERBOSE>=2) std::cout<<"Initial vertex position = "<<initialVertex[icand][0]<<","<<initialVertex[icand][1]<<","<<initialVertex[icand][2]<<","<<initialVertex[icand][3]<<", Candidate vertex position = ("<<tCandRecoVtxPos[0]<<","<<tCandRecoVtxPos[1]<<","<<tCandRecoVtxPos[2]<<","<<tCandRecoVtxPos[3]<<", NLL = "<<fmin<<std::endl;
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
	}
	
	delete minimizer;

	tBestReconstructedVertexPosition = tReconstructedVertexPosition[0]; 
		
	clock_t timeEnd=clock();
	
	if(verbose) 
		std::cout<<"MinimizeVertex: Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<std::endl<<std::endl;

	return tReconstructedVertexPosition;
	//return tBestReconstructedVertexPosition;
}


std::vector< std::vector<double> > LEAF::MinimizeVertex_Main(	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool likelihood,bool average,double lowerLimit, double upperLimit, int directionality){
		
		
	std::vector<std::thread> lThreadList;
	std::vector< std::vector<double> > lOutputFinal;
	
	int iCand_Step = nCandidates / fThread;
	
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
void LEAF::MinimizeVertex_thread(	
		int iStart, int iIte,	
		std::vector< std::vector<double> > initialVertex, double * limits, double stepSize, int nhits,
		int nCandidates, int tolerance, int verbose,bool /*likelihood*/,bool /*average*/,double lowerLimit, double upperLimit, int directionality){
		
	int iEnd = (iIte+iStart);
	if ( iEnd > nCandidates ) iEnd = nCandidates;
	
	// Set local TGraph2D

	std::vector<struct FitPosition>  tVtxContainer;

	mtx.lock();
	TFitter * minimizer = new TFitter(8);//4=nb de params?
	mtx.unlock();
	TMinuit * minuit = minimizer->GetMinuit();
	
	double arglist[20];
	int err=0;
	//arglist[0]=0;
	double p1=verbose-1;
	minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
	minuit->SetErrorDef(1);
	//minuit->SetPrintLevel(verbose-1); // quiet mode
	if(VERBOSE<2){
		minuit->mnexcm("SET NOWarnings",0,0,err);
	}
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
		if(VERBOSE>=2) std::cout<<"Initial vertex position = "<<initialVertex[icand][0]<<","<<initialVertex[icand][1]<<","<<initialVertex[icand][2]<<","<<initialVertex[icand][3]<<", Candidate vertex position = ("<<tCandRecoVtxPos[0]<<","<<tCandRecoVtxPos[1]<<","<<tCandRecoVtxPos[2]<<","<<tCandRecoVtxPos[3]<<", NLL = "<<fmin<<std::endl;
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
	delete minimizer;
	mtx.unlock();
}

struct LEAF::FitterOutput LEAF::MakeFit(const HitCollection<Hit>* lHitCol, const TimeDelta lTriggerTime, bool bMultiPMT) {

	fHitCollection  = lHitCol; 
	fTriggerTime    = lTriggerTime;
	fTimeCorrection = fHitCollection->timestamp - fTriggerTime;

	// Get Hits total

	//int iHitsTotal = fHitInfo.size();
	int iHitsTotal = fHitCollection->Size();
	
	if ( fPositionList.size() == 0 ) {
		// Make position lists
      		this->MakePositionList();
	}
	
	bool bDirectionnality = false;
	if ( bMultiPMT ) {
		bDirectionnality = fUseDirectionality;
	}
	
	//Proceed to the fit.
	double * tLimits = new double[4];
	
	struct FitterOutput fOutput;
	
	fOutput.Vtx[0]		= 0.;
	fOutput.Vtx[1]		= 0.;
	fOutput.Vtx[2]		= 0.;
	fOutput.Vtx[3]		= 0.;
	fOutput.NLL		= 0.;
	
	fOutput.InTime		= 0;
	
	fOutput.True_NLLDiff	= 0;
	fOutput.True_TimeDiff	= 0;
	fOutput.True_TistDiff	= 0;
	
	fLastLowerLimit       = 0.;
	fLastUpperLimit       = 0.;
		 	
	if ( true ) {
	
		if(iHitsTotal>0){
		
			TStopwatch timer;
			timer.Reset();
			timer.Start();
				
			//2.a. Coarse vertex search from scratch
			double dStepSize=fSearchVtxStep;//in centimeters
			//std::vector< std::vector<double> > tRecoVtxPos = this->SearchVertex(iHitsTotal,fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
			std::vector< std::vector<double> > tRecoVtxPos = this->SearchVertex_Main(iHitsTotal,fSearchVtxTolerance,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,false);
			
						
			if(VERBOSE>=2){
				for(unsigned int a=0;a<tRecoVtxPos.size();a++) std::cout<<"Candidate vertex [ " << a << " ] = " <<tRecoVtxPos[a][0] << " , " << tRecoVtxPos[a][1] << " , " << tRecoVtxPos[a][2] << " , " << tRecoVtxPos[a][3] << ", NLL = "<< tRecoVtxPos[a][4] << std::endl;
			}
			timer.Stop();
			//std::cout << "SearchVertex took: " << timer.RealTime() << std::endl;
			
			////////////////////////////////////////
					
			double dStepSizeFinal=10;
			int iToleranceFinal=1;
			if(fStepByStep){
				//2.b. Refined vertex search around candidate found previously
				// -> Commented since ages, removed on 2020/11/20
				
				////////////////////////////////////////
		    
				//2.c. More refined vertex search around candidate found previously
				double dStepSizeFine2=100;
				int iToleranceFine2=20;
				for(int i=0;i<4;i++) tLimits[i] = dStepSize;//particleStart[i]*1e-2;//onversion to meter
				std::vector< std::vector<double> > tRecoVtxPosFine2 = this->SearchVertexFine(tRecoVtxPos,tLimits,dStepSizeFine2,iHitsTotal,fSearchVtxTolerance,iToleranceFine2,VERBOSE,true,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,bDirectionnality);
				//////////////////////////////////////////
		    
				//2.d. Final refined vertex search around candidate found previously
				double dStepSizeFine3=20;
				int iToleranceFine3=1;
				for(int i=0;i<4;i++) tLimits[i] = dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
				std::vector< std::vector<double> > tRecoVtxPosFine3 = this->SearchVertexFine(tRecoVtxPosFine2,tLimits,dStepSizeFine3,iHitsTotal,iToleranceFine2,iToleranceFine3,VERBOSE,true,false,fSTimePDFLimitsQueueNegative,fSTimePDFLimitsQueuePositive,bDirectionnality);
				if(VERBOSE == 1){
					for(int a=0; a<iToleranceFine3; a++) std::cout<<"Final candidate vertex time = "<<tRecoVtxPosFine3[a][0]<<", x = "<<tRecoVtxPosFine3[a][1]<<", y = "<<tRecoVtxPosFine3[a][2]<<", z="<<tRecoVtxPosFine3[a][3]<<std::endl;
				}
				//////////////////////////////////////////	
				
				for(int i=0;i<4;i++) tLimits[i] = 100;//dStepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
				fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPosFine3,tLimits,dStepSizeFinal,iHitsTotal,iToleranceFine3,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,bDirectionnality);
				if(VERBOSE >= 1){
					for(int a=0;a<iToleranceFinal;a++){
						std::cout<<"Final candidate vertex time = "<<fRecoVtxPosFinal[a][0]<<", x = "<<fRecoVtxPosFinal[a][1]<<", y = "<<fRecoVtxPosFinal[a][2]<<", z="<<fRecoVtxPosFinal[a][3]<<std::endl;
					}
				}
			}
			else {
			
				timer.Reset();
				timer.Start();
				
				for(int i=0;i<4;i++) tLimits[i] = 2*dStepSize;
				
				//fRecoVtxPosFinal = this->MinimizeVertex(tRecoVtxPos,tLimits,dStepSizeFinal,iHitsTotal,fSearchVtxTolerance,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,bDirectionnality);
				fRecoVtxPosFinal = this->MinimizeVertex_Main(tRecoVtxPos,tLimits,dStepSizeFinal,iHitsTotal,fSearchVtxTolerance,iToleranceFinal,VERBOSE,true,false,fMinimizeLimitsNegative,fMinimizeLimitsPositive,bDirectionnality);
				
				timer.Stop();
				//std::cout << "Minimizer took: " << timer.RealTime() << std::endl;
			}
			
			//3.a. Test the NLL of the true vertex to compare with reco
			// -> Commented since ages, removed on 2020/11/20

			
			int iInTime = 0;
			
			for(int ihit = 0; ihit < iHitsTotal; ihit++){
				//Hit lHit = fHitInfo[ihit];
				Hit lHit = fHitCollection->At(ihit);
				
				int iPMT = lHit.PMT;
				PMTInfo lPMTInfo = (*fPMTList)[iPMT];
				
				double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
				double distance = Astro_GetDistance(lPMTInfo.Position,fRecoVtxPosFinal[0]);
				double tof = distance / fLightSpeed;
				double residual = hitTime - tof - fRecoVtxPosFinal[0][3];
				if(residual > fSTimePDFLimitsQueueNegative && residual < fSTimePDFLimitsQueuePositive){
					iInTime++;
				}
			}
			
			// Fill Output
			fOutput.Vtx[0] 	= fRecoVtxPosFinal[0][0];
			fOutput.Vtx[1] 	= fRecoVtxPosFinal[0][1];
			fOutput.Vtx[2] 	= fRecoVtxPosFinal[0][2];
			fOutput.Vtx[3]   	= (fRecoVtxPosFinal[0][3] - fTimeCorrection) / TimeDelta::ns;
			fOutput.NLL 		= fRecoVtxPosFinal[0][4];  
			
			fOutput.InTime 	= iInTime;
			
			//std::cout << " Final vertex " << fOutput.Vtx[0] << " " << fOutput.Vtx[1] << " " << fOutput.Vtx[2] << " Time " << fOutput.Vtx[3] << std::endl;
		}
	}
	
	return fOutput;
}
