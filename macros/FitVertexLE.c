//Documentation by B.Quilain, 2019/10.23
//V1 was the first vertex fitter, adapted to WCSim output in Spring 2019, i.e. where 2 PMT types digitizer results were stored into a different event number i.e. one physics event #0 was digitized into 2 digitizers (1/PMT type) and each results stored into event#0 (PMT type 0) and event#1 (PMT type 1).
//V2 is required because B.Quilain changed WCSim output to store the 2 digitizer results into the same event, but in a different branch name:
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TBranch.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TFitter.h>
#include <TMinuit.h>
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
WCSimRootGeom *geo = 0; 
#define KILLFAKE
//#define KILLHALF
TRandom3 * trand = new TRandom3();
TRandom3 * trand2 = new TRandom3();

bool maskPMT=false;//true;//false;//true;
int pmtTypeMasked=1;
double maskedProbability=0.4;

bool isHE=false;//true;
bool vertexFit=true;
bool stepByStep=false;
bool useDirectionality=true;//false;//true;//true;//false;//true;//false;//true;//false;//true;//true;//true;//true;//false;//true;//false;//true;//false;//true;//false;//true;//false;
bool limitmPMT = true;//false;//true;//false;//true;//true;//false;//true;//If true, we apply same cuts for timing as B&L PMT when searching vertex. Otherwise, we do not apply them and use all the hits. The latter is particularly useful when directionality is added, as it can kill DR hits.
int nAveraging=20;//Number of points we used to average the vertex on.
double timeWindowSizeFull=1500;
//#define

double ** essentialHitInfo;
double ** essentialHitInfo_pmtType0;//Kepp a copy of the first PMT type list of hits
int nhits_pmtType0;

float cvacuum = 3e8*1e2 / 1e9;//speed of light, in centimeter per ns.
float nindex = 1.373;//1.385;//1.373;//refraction index of water
double lightSpeed = cvacuum/nindex;

using namespace std;
// Simple example of reading a generated Root file
const int nPMTtypes = 2;
const int npmts_in_mpmt = 19;
const int nGroupsPMTs = 3;//Different PMT config in an mPMT. I assumed here a rotational symetry of the mPMT, so PMT 1 to 12 are the same, 13 to 18 are the same and 19 is separated
const int maxPMTTubes = 1e6;

//If we use the Bayesian CL for directionality
double DR_dir_proba[nPMTtypes][nGroupsPMTs] ={{0.}};

double DRns[nPMTtypes]={8400*2e4*1e-9,100*19*5e3*1e-9};
double ftimePDFLimits=100;
double stimePDFLimits=8;
double hitTimeLimitsNegative=-5;//-10;//-5;//-5
double hitTimeLimitsPositive=7;//15;//7;//7
double stimePDFLimitsQueueNegative=-3;//-5;//-3;//-1.5;//-3;//-10;
double stimePDFLimitsQueuePositive=4;//10;//4;//3;//4;//500;
double stimePDFLimitsQueueNegative_fullTimeWindow=0;//-1.5;//-3;//-10;
double stimePDFLimitsQueuePositive_fullTimeWindow=0;//3;//4;//500;
double PDFnormalization_fullTimeWindow=0;//likelihood?splineIntegral(stimePDFQueue[pmtType],stimePDFLimitsQueueNegative_fullTimeWindow,stimePDFLimitsQueuePositive_fullTimeWindow):1;

//2D plots
TH2D * ChargeProfile2D[nPMTtypes];
TH2D * ChargeProfile2D_onlyFront[nPMTtypes];
TH2D * ChargeProfile2DTop[nPMTtypes];
TH2D * ChargeProfile2DCylindricCoordinate[nPMTtypes];
TH2D * ChargeProfile2DRelativeAngle[nPMTtypes];
TH2D * TimeTOFProfileXTOF[nPMTtypes];

//1D plots of charge & hit
TH1D * ChargeProfile[nPMTtypes];
TH1D * HitProfile[nPMTtypes];

//1D time profiles
TH1D * TimeProfile[nPMTtypes];
TH1D * TimeHitProfile[nPMTtypes];
TH1D * TimeTOFProfile[nPMTtypes];
TH1D * HitTimeTOFProfile[nPMTtypes];

//Charge per PMT
TH1D * ChargePerPMT[nPMTtypes];

//Total charge
TH1D * TotalCharge[nPMTtypes];
TH1D * TotalHit[nPMTtypes];
TH2D * TotalChargeXdWall[nPMTtypes];

TH3D * timeWindowXsignal[nPMTtypes];
TH3D * timeWindowXDR[nPMTtypes];
TH1D * TimeTOFProfile_currentEvent[nPMTtypes];

//TH1D * EfficiencyLE[nPMTtypes];
TH1D * fitTimeToTrue[nPMTtypes];
TH1D * fitDistanceToTrue[nPMTtypes];
TH2D * fit4DdistanceXNLLdifference[nPMTtypes];
TH2D * fitNLLprofile_currentEvent[nPMTtypes];
TH2D * normfitNLLprofile_currentEvent[nPMTtypes];

TSpline3 * stimePDF[nPMTtypes];TF1 * ftimePDF[nPMTtypes];
TSpline3 * stimePDFQueue[nPMTtypes];
TSpline3 * stimePDFDR[nPMTtypes];
//TH2D * hPMTDirectionality[nPMTtypes][npmts_in_mpmt];
TH1D * hPMTDirectionality_1D[nPMTtypes][nGroupsPMTs];


double * crossProduct(double* v0, double* v1, double * v2){
  v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
  v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
  v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
  return v2;
}

double groupPMTs(int pmtType, int pmt_number){
  if(pmtType==0) return 0;
  else{
    if(pmt_number<=12) return 0;
    else if(pmt_number<=18) return 1;
    else return nGroupsPMTs-1;
  }
}

void vectorVertexPMT(double * vertex, int tubeNumber, int pmtType, double * vDirection,int pmt_number_in_mpmt,int verbose=0){
  clock_t timeStart=clock();
  const int limit_maxPMT_per_mPMT = 50;
  double pmtCenter[3]={0.};
  double pmt_position_local[3]={0.};
  int nPMT_per_mPMT=0;
  double PMTpos[3];
  double PMTdir[3];
  WCSimRootPMT pmt;
  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
  else pmt  = geo->GetPMT(tubeNumber-1,true);
  for(int j=0;j<3;j++){
    PMTpos[j] = pmt.GetPosition(j);
    PMTdir[j] = pmt.GetOrientation(j);
  }
  int mPMT_number = pmt.GetmPMTNo();
  pmt_number_in_mpmt = pmt.GetmPMT_PMTNo();
  double particleRelativePMTpos[3];
  for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vertex[j];
  if(verbose == 3){
    cout<<"Original PMT number = "<<pmt_number_in_mpmt<<endl;  
    cout<<"Direction = "<<PMTdir[0]<<", "<<PMTdir[1]<<", "<<PMTdir[2]<<endl;
  }
  for(int t = -limit_maxPMT_per_mPMT/2;t<limit_maxPMT_per_mPMT/2;t++){
    int tubeNumber2=tubeNumber+t;
    if(tubeNumber2<1) continue;
    
    WCSimRootPMT pmt2;
    if(pmtType == 0) pmt2 = geo->GetPMT(tubeNumber2-1,false);
    else pmt2  = geo->GetPMT(tubeNumber2-1,true);
    int mPMT_number2 = pmt2.GetmPMTNo();
    int pmt_number_in_mpmt2 = pmt2.GetmPMT_PMTNo();
    
    if(mPMT_number2 == mPMT_number){
      nPMT_per_mPMT++;
      if(verbose == 3) cout<<"PMT number = "<<pmt_number_in_mpmt2<<", ";  
      for(int j=0;j<3;j++){
	pmtCenter[j]+=pmt2.GetPosition(j);
	if(verbose == 3) cout<<pmt2.GetPosition(j)<<", ";
      }
      if(verbose == 3) cout<<endl;
    }
  }
  
  if(verbose == 3) cout<<"Position with respect to mean = ";
  for(int j=0;j<3;j++){
    if(nPMT_per_mPMT!=0) pmtCenter[j]/=nPMT_per_mPMT;
    pmt_position_local[j]= PMTpos[j] - pmtCenter[j];
    if(verbose == 3) cout<<pmt_position_local[j]<<", ";
  }
  if(verbose == 3) cout<<endl;
  
  //Now, let's check the angle from the PMT to the center in the referential orthogonal to the PMT direction.
  //PMT direction is given by PMTdir
  //We can construct the vector orthogonal to PMTdir and the direction from PMT to center. We will then normalize it.
  //And finally, we should use that vertex and PMT dir to construct the orthogonal vector. That last vector will be in the same plane than PMTdir and PMT to center.
  double * vCenter = new double[3];
  double norm=0;
  for(int j=0;j<3;j++){
    vCenter[j]=-pmt_position_local[j];//Minus since we want the vector from PMT to center (and not opposite).
    norm+=pow(vCenter[j],2);
  }
  for(int j=0;j<3;j++){
    if(norm!=0) vCenter[j]/=TMath::Sqrt(norm);
  }
  
  double ** v = new double*[3];
  for(int x=0;x<3;x++) v[x] = new double[3];
  norm=0;
  for(int j=0;j<3;j++){
    v[2][j]=PMTdir[j];
    norm+=pow(v[2][j],2);
  }
  for(int j=0;j<3;j++){
    if(norm!=0) v[2][j]/=TMath::Sqrt(norm);
  }
  
  norm=0;
  crossProduct(v[2],vCenter,v[1]);
  for(int j=0;j<3;j++){
    norm+=pow(v[1][j],2);
  }
  for(int j=0;j<3;j++){
    if(norm!=0) v[1][j]/=TMath::Sqrt(norm);
  }
  
  norm=0;
  crossProduct(v[1],v[2],v[0]);
  for(int j=0;j<3;j++){
    norm+=pow(v[0][j],2);
  }
  for(int j=0;j<3;j++){
    if(norm!=0) v[0][j]/=TMath::Sqrt(norm);
  }

  //
  //Then, we should project the PMT to vertex in this new basis
  double * vPMTvertex = new double[3];
  for(int x=0;x<3;x++){
    double ps=0;//scalar product
    for(int j =0;j<3;j++){
      ps+=(-particleRelativePMTpos[j])*v[x][j];
      //ps+=(-PMTpos[j])*v[x][j];
    }
    vPMTvertex[x]=ps;
  }
	  
  if(verbose == 3){
    cout<<"v dir:"<<endl;
    for(int x=0;x<3;x++){
      cout<<"vdir["<<x<<"]: ";
      cout<<v[x][0]<<", "<<v[x][1]<<", "<<v[x][2]<<endl;
    }
    cout<<"To center PMT = "<<vCenter[0]<<", "<<vCenter[1]<<", "<<vCenter[2]<<endl;
    cout<<"PMT vertex = "<<vPMTvertex[0]<<", "<<vPMTvertex[1]<<", "<<vPMTvertex[2]<<endl;
  }
  //First, get radius of this vector:
  double radius=0;
  for(int j =0;j<3;j++){
    radius+=pow(vPMTvertex[j],2);
  }
  radius=TMath::Sqrt(radius);
  //We just want phi:
  for(int j =0;j<3;j++){
    vPMTvertex[j]/=radius;
  }
  double cosTheta = vPMTvertex[2];
  double tanPhi = vPMTvertex[1]/vPMTvertex[0];
  //atan is defined between -pi/2 and pi/2
  double Theta = TMath::ACos(cosTheta);
  double Phi;
  if(vPMTvertex[0]>=0) Phi=TMath::ATan(tanPhi);
  else Phi=-TMath::ATan(tanPhi);
  //The phi angle is the angle of
  //We can first construct the cosine
  vDirection[0]=Phi*180/TMath::Pi();
  vDirection[1]=Theta*180/TMath::Pi();

  if(verbose == 3){
    cout<<"Number of PMTs in mPMT="<<nPMT_per_mPMT<<", PMT #"<<pmt_number_in_mpmt<<", position local = "<< pmt_position_local[0] << ", " << pmt_position_local[1] << ", " << pmt_position_local[2] << endl;
    cout<<"Phi angle="<<vDirection[0]<<", theta = "<<vDirection[1]<<endl;
  }
  
  for(int x=0;x<3;x++) delete v[x];
  delete vCenter;
  delete v;
  delete vPMTvertex;
  clock_t timeEnd=clock();
  if(verbose == 3) cout<<"Time for the loop = "<<timeEnd-timeStart<<endl;
  //delete vPMTvertex;
  
}

double findDirectionTheta(double * vertex,int tubeNumber,int pmtType,int pmt_number_in_mpmt,int verbose=0){
  clock_t timeStart=clock();
  double PMTpos[3];
  double PMTdir[3];
  WCSimRootPMT pmt;
  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
  else pmt  = geo->GetPMT(tubeNumber-1,true);
  for(int j=0;j<3;j++){
    PMTpos[j] = pmt.GetPosition(j);
    PMTdir[j] = pmt.GetOrientation(j);
  }

  int mPMT_number = pmt.GetmPMTNo();
  pmt_number_in_mpmt = pmt.GetmPMT_PMTNo();
  double particleRelativePMTpos[3];
  for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - vertex[j];
  if(verbose == 3){
    cout<<"Original PMT number = "<<pmt_number_in_mpmt<<endl;  
    cout<<"Direction = "<<PMTdir[0]<<", "<<PMTdir[1]<<", "<<PMTdir[2]<<endl;
  }

  double norm=0;
  double ps=0;
  for(int j =0;j<3;j++){
    ps+=(-particleRelativePMTpos[j])*PMTdir[j];
    norm+=pow(particleRelativePMTpos[j],2);
  }
  double cosTheta=ps/TMath::Sqrt(norm);
  double Theta=TMath::ACos(cosTheta)*180/TMath::Pi();
  
  if(verbose == 3){
    cout<<"PMT #"<<pmt_number_in_mpmt<< ", Theta = "<< Theta << endl;
  }
  
  clock_t timeEnd=clock();
  if(verbose == 3) cout<<"Time for the loop = "<<timeEnd-timeStart<<endl;
  return Theta;
  //delete vPMTvertex;
  
}

double findNLLDirectionality(double vertexPosition[4],int nhits,int verbose,double lowerLimit,double upperLimit){
  double NLL = 0;
  double *info;
  double PDFnormalization=1;//splineIntegral(stimePDFQueue[pmtType],lowerLimit,upperLimit);
  double * vPos = new double[3];
  double * vDirection = new double[2];//Return phi and theta.

  for(int ihit = 0; ihit < nhits; ihit++){
    info = essentialHitInfo[ihit];
    double hitTime = info[0];
    double hitPosition[3];
    double distance = 0;
    int pmtType=info[4];
    int tubeNumber=info[5];
    int pmt_number_in_mpmt;

    for(int i=0;i<3;i++){
      vPos[i] = vertexPosition[i+1];
    }

    if(pmtType==0) continue;
    
    for(int i=0;i<3;i++){
      hitPosition[i] = info[i+1] - vertexPosition[i+1];
      distance+=pow(hitPosition[i],2);
    }
    distance=TMath::Sqrt(distance);
    double tof = distance / lightSpeed;
    double residual = hitTime - tof - vertexPosition[0];
    bool condition;
    if(pmtType == 0 || (pmtType == 1 && limitmPMT) ) condition = residual > lowerLimit && residual < upperLimit;
    else condition = true;
    if(condition){
      //if(1/*residual > lowerLimit && residual < upperLimit*/){//Do not use it as fluctuating the number of hits used in one event does allow a fair comparison with this likelihood. If we want to restrict the likelihood, use the "DirectionalNLLBayes".
      //vectorVertexPMT(vPos,tubeNumber,pmtType,vDirection,pmt_number_in_mpmt,verbose);
      //proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(vDirection[1]));
      double theta = findDirectionTheta(vPos,tubeNumber,pmtType,pmt_number_in_mpmt,verbose);
      int pmtGroup=groupPMTs(pmtType,pmt_number_in_mpmt);
      double proba;
      proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(theta));

      
      if(proba==0){
	double min=0;
	for(int ibinx=hPMTDirectionality_1D[pmtType][pmtGroup]->GetNbinsX();ibinx>=1;ibinx--){//Search the last bin of the distri non-zero and use its value
	  if(hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx) !=0 ){
	    min = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx);
	    break;
	  }
	}
	//proba=min;
	proba=1e-20;
	cout<<"We are at proba = 0, theta = "<<theta<<", proba used = "<<proba<<endl;
      }
      NLL += -TMath::Log(proba);
      
      if(verbose==3){
	//cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << vDirection[1] << ", proba="<<proba<<endl;
	cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << theta << ", proba="<<proba<<endl;
      }
    }
  }
  if(verbose==2) cout<<"NLL directionnel="<<NLL<<endl;
  delete vDirection, vPos;
  return NLL;
}

//Created to compare 2 hypotheses: DR and signal.
//Useful when we have a number of hits that fluctuates between 2 candidate vertex. The main example is when we restrict the likelihood to a given time window
//Indeed, in that case, as the PDF <=1, the more hits we have, the less the Likelihood. But, what we want is to be able to compare cases with other vertices which corresponds to different number of hits in time window.
double findNLLDirectionalityBayes(double vertexPosition[4],int nhits,int verbose,double lowerLimit,double upperLimit){
  double NLLsignal = 0.9;
  double NLLdr = 0.1;
  if(!limitmPMT){//Proba of having DR in a large time window is higher
    NLLsignal = 1e-6;
    NLLdr = 1-1e-6;
  }
  double *info;
  double PDFnormalization=1;//splineIntegral(stimePDFQueue[pmtType],lowerLimit,upperLimit);
  double * vPos = new double[3];
  double * vDirection = new double[2];//Return phi and theta.

  for(int ihit = 0; ihit < nhits; ihit++){
    info = essentialHitInfo[ihit];
    double hitTime = info[0];
    double hitPosition[3];
    double distance = 0;
    int pmtType=info[4];
    int tubeNumber=info[5];
    int pmt_number_in_mpmt;

    for(int i=0;i<3;i++){
      vPos[i] = vertexPosition[i+1];
    }

    if(pmtType==0) continue;


    for(int i=0;i<3;i++){
      hitPosition[i] = info[i+1] - vertexPosition[i+1];
      distance+=pow(hitPosition[i],2);
    }
    distance=TMath::Sqrt(distance);
    double tof = distance / lightSpeed;
    double residual = hitTime - tof - vertexPosition[0];
    //if(1){
    bool condition;
    if(pmtType == 0 || (pmtType == 1 && limitmPMT) ) condition = residual > lowerLimit && residual < upperLimit;
    else condition = true;
    if(condition){
      //vectorVertexPMT(vPos,tubeNumber,pmtType,vDirection,pmt_number_in_mpmt,verbose);
      //proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(vDirection[1]));
      double theta = findDirectionTheta(vPos,tubeNumber,pmtType,pmt_number_in_mpmt,verbose);
      int pmtGroup=groupPMTs(pmtType,pmt_number_in_mpmt);
      double proba;
      proba = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(hPMTDirectionality_1D[pmtType][pmtGroup]->FindBin(theta));
      
      if(proba==0){
	double min=0;
	for(int ibinx=hPMTDirectionality_1D[pmtType][pmtGroup]->GetNbinsX();ibinx>=1;ibinx--){//Search the last bin of the distri non-zero and use its value
	  if(hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx) !=0 ){
	    min = hPMTDirectionality_1D[pmtType][pmtGroup]->GetBinContent(ibinx);
	    break;
	  }
	}
	proba=min;
	cout<<"We are at proba = 0, theta = "<<theta<<", proba used = "<<proba<<endl;
      }
      //NLL += -TMath::Log(proba);
      NLLsignal*=proba;
      NLLdr*=DR_dir_proba[pmtType][pmtGroup];

      if(verbose==3){
	//cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << vDirection[1] << ", proba="<<proba<<endl;
	cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << theta << ", proba="<<proba<<endl;
      }
    }
    else if(verbose == 2) cout<<"Out of boundaries for directional L"<<endl;
  }
  double NLL = NLLsignal / (NLLsignal+NLLdr);
  if(NLL>0) NLL=-TMath::Log(NLL);
  if(verbose==2) cout<<"NLL directionnel Bayesian="<<NLL<<endl;
  delete vDirection, vPos;
  return NLL;
}

Double_t valueTimePDF(Double_t *time, Double_t *par){
  if(time[0]<par[0] || time[0]>par[1]) return 0;
  else{
    int type=((int) par[2]);
    double value = stimePDF[type]->Eval(time[0]);
    if(value<0) return 0;//May arrive for few points which are e.g. at -1e-25, due to imperfect extrapolation
    else stimePDF[type]->Eval(time[0]);
  }
}

void loadSplines(){
  TFile *fSplines, *fSplines2;
  if(isHE){
    fSplines = new TFile("../inputs/timePDF_HE.root","read");//To generate with my code ProduceWSPlots.c
    fSplines2 = fSplines;
  }
  else{
    fSplines = new TFile("../inputs/timePDF_DRnew_Large.root","read");//To generate with my code ProduceWSPlots.c
    fSplines2 = new TFile("../inputs/timePDF_Directionality.root","read");//To generate with my code ProduceWSPlots.c
  }
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    //Load 1D t-tof splines
    stimePDFQueue[pmtType] = (TSpline3*) fSplines->Get(Form("splineExpoQueue%d_%d",0,pmtType));
    stimePDFDR[pmtType] = (TSpline3*) fSplines->Get(Form("splineDR%d_%d",0,pmtType));

    stimePDFLimitsQueueNegative_fullTimeWindow = stimePDFQueue[pmtType]->GetXmin();
    stimePDFLimitsQueuePositive_fullTimeWindow = stimePDFQueue[pmtType]->GetXmax();
    cout<<"Min="<<stimePDFLimitsQueueNegative_fullTimeWindow<<", max="<<stimePDFLimitsQueuePositive_fullTimeWindow<<endl;
    if(stimePDFLimitsQueueNegative<stimePDFLimitsQueueNegative_fullTimeWindow) stimePDFLimitsQueueNegative=stimePDFQueue[pmtType]->GetXmin();//To avoid to use PDF where it is not defined
    if(stimePDFLimitsQueuePositive>stimePDFLimitsQueuePositive_fullTimeWindow) stimePDFLimitsQueuePositive=stimePDFQueue[pmtType]->GetXmax();//same here.

    //Load 3D directionality histograms.
    for(int pmtGroup=0;pmtGroup<nGroupsPMTs;pmtGroup++){
      hPMTDirectionality_1D[pmtType][pmtGroup] = (TH1D*) fSplines2->Get(Form("hPMTDirectionality_1D_%d_%d_%d",0,pmtType,pmtGroup));
      DR_dir_proba[pmtType][pmtGroup] = hPMTDirectionality_1D[pmtType][pmtGroup]->Integral() / hPMTDirectionality_1D[pmtType][pmtGroup]->GetNbinsX();
      cout<<hPMTDirectionality_1D[pmtType][pmtGroup]->GetMean()<<", and DR = "<<DR_dir_proba[pmtType][pmtGroup]<<endl;
    }
    /*
    bohPMTDirectionality_1Dol testSplines=true;
    if(testSplines){
      cout<<"testing spline for pmt type = " << pmtType << endl;
      cout<<"Min value of the spline range = " << stimePDF[pmtType]->GetXmin() << endl;
      for(double t=-stimePDFLimits;t<stimePDFLimits;t+=1){
	cout << "Time " << t << "ns, value = "<< stimePDF[pmtType]->Eval(t) << endl;
	//cout << "Time " << t << "ns, value = "<< ftimePDF[pmtType]->Eval(t) << endl;
      }
      }
    ftimePDF[pmtType] = new TF1(Form("ftimePDF%d",pmtType),valueTimePDF,-ftimePDFLimits,ftimePDFLimits,3);
    ftimePDF[pmtType]->SetParameter(0,-20);
    ftimePDF[pmtType]->SetParameter(1,20);
    ftimePDF[pmtType]->SetParameter(2,pmtType);*/
  }
}
double splineIntegral(TSpline3 * s,double start,double end,double stepSize=5e-1){
  double integral=0;
  for(double i=start;i<end;i+=stepSize){
    integral+=stepSize*s->Eval(i+stepSize/2);
  }
  return integral;
}
double splineIntegralAndSubstract(TSpline3 * s0,TSpline3 * s1,double start,double end,double stepSize=5e-1){
  double integral=0;
  for(double i=start;i<end;i+=stepSize){
    integral+=stepSize*(s0->Eval(i+stepSize/2)-s1->Eval(i+stepSize/2));
  }
  return integral;
}

double splineIntegralExpo(TSpline3 * s,double start,double end,double sigma,double stepSize=5e-1){
  double integral=0;
  for(double i=start;i<end;i+=stepSize){
    integral+=stepSize*s->Eval(i+stepSize/2)*TMath::Gaus(i+stepSize/2,0,sigma);
  }
  return integral;
}

double findNLL(double vertexPosition[4],int nhits,bool likelihood,int verbose,double lowerLimit, double upperLimit, bool killEdges=false, bool scaleDR=false,int directionality=false){
  double NLL = 0;
  double *info;
  double PDFnormalization=1;//splineIntegral(stimePDFQueue[pmtType],lowerLimit,upperLimit);
  //if(PDFnormalization_fullTimeWindow!=0) PDFnormalization/=PDFnormalization_fullTimeWindow;
  //We know the DR integral over a period of time. But the proba to have DR in one bin depends on the binning of the Spline, since it depends on the bin size.
  //Shouldn't we: 1. Extract the graph from the spline3
  /* TH1D* h = (TH1D*) ->GetHistogram();
  double signalPEinTime=;
  double DRinTime=DRns[pmtType]*(upperLimit-lowerLimit);
  //We should adapt signal / DR the same way as in time window total
  double signalPE = std::max(nhits - timeWindowSizeFull*DRns[pmtType],0.);
  double signalIncrease=signalPE/signalPEInTime;//We will scale the signal to have the same integral as number of hits.
  //In that situation, it is simple to deduce the DR value: in a time window 
  double DRtotal = timeWindowSizeFull*DRns[pmtType];
  */
  //double DRtotal = DRns[pmtType]
  //double DR=;//To rescale the PDF, we should know the relative integral of signal vs DR
  double nhits_perType[nPMTtypes];
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    nhits_perType[pmtType]=0;
  }
  for(int ihit = 0; ihit < nhits; ihit++){
    double *info_temp = essentialHitInfo[ihit];
    int pmtType=info_temp[4];
    nhits_perType[pmtType]++;
  }
  double signalOverNoise[nPMTtypes],DRIntegral[nPMTtypes],signalIntegral[nPMTtypes];
  if(scaleDR){
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      double signalDR, signalPE;
      signalDR=timeWindowSizeFull*DRns[pmtType];//Over the whole time window
      signalPE=std::max(nhits_perType[pmtType] - signalDR,0.);//Over the whole time window.
      double signalPETime=signalPE;//I assume that all signal is here.
      double signalDRTime=(upperLimit-lowerLimit)*signalDR/timeWindowSizeFull;
      signalOverNoise[pmtType]=signalPETime/signalDRTime;//Over the 
      DRIntegral[pmtType]=splineIntegral(stimePDFDR[pmtType],lowerLimit,upperLimit);
      signalIntegral[pmtType]=splineIntegralAndSubstract(stimePDFQueue[pmtType],stimePDFDR[pmtType],lowerLimit,upperLimit);
      if(verbose==3) cout<<"nhits="<<nhits_perType[pmtType]<<", DR average="<<signalDR<<", signal over noise="<<signalOverNoise[pmtType]<<", signal integral="<<signalIntegral[pmtType]<<", DR integral="<<DRIntegral[pmtType]<<",in integral="<<signalIntegral[pmtType]/DRIntegral[pmtType]<<endl;
    }
    /*    if(verbose==2){
      cout<<"Lower limit = "<<lowerLimit<<", proba = "<<stimePDFQueue[pmtType]->Eval(lowerLimit)<<endl;
      cout<<"Upper limit = "<<upperLimit<<", proba = "<<stimePDFQueue[pmtType]->Eval(upperLimit)<<endl;
      }*/
  }

  
  for(int ihit = 0; ihit < nhits; ihit++){
    info = essentialHitInfo[ihit];
    double hitTime = info[0];
    double hitPosition[3];
    double distance = 0;
    int pmtType=info[4];
    for(int i=0;i<3;i++){
      hitPosition[i] = info[i+1] - vertexPosition[i+1];
      distance+=pow(hitPosition[i],2);
    }
    distance=TMath::Sqrt(distance);
    double tof = distance / lightSpeed;
    double residual = hitTime - tof - vertexPosition[0];
    double proba;
#ifdef KILLHALF
    if(pmtType==1){
      double t=trand->Uniform(0,1);
      if(t>0.5) continue;
    }
#endif

    if(likelihood){
      bool condition;
      //if(pmtType == 0 || (pmtType == 1 && limitmPMT) ) condition = residual > lowerLimit && residual < upperLimit;
      //else condition = true;
      condition = residual > lowerLimit && residual < upperLimit;
      if(condition){
	//if(pmtType==1) continue;
	proba = stimePDFQueue[pmtType]->Eval(residual);
	//else proba = 1;//stimePDFQueue[pmtType]->Eval(residual);
	if(verbose==3){
	  cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[0] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<endl;
	  cout<<"Residual="<<residual<<", proba="<<proba<<", pmt type="<<pmtType<<endl;
	}
	if(scaleDR){
	  //First, substract the DR from the PDF to keep only the signal:
	  double DR=stimePDFDR[pmtType]->Eval(residual);
	  proba-=DR;
	  //Then, scale signal so that signal integral / DR integral = signalOverNoise
	  //To do so, we should scale signal so that integral = DR integral * signalOverNoise
	  //And of course, to rescale, we should divide by the signal integral
	  double Factor = DRIntegral[pmtType] * signalOverNoise[pmtType] / signalIntegral[pmtType];
	  proba*=Factor;
	  //And add again the DR
	  proba+=DR;
	  if(verbose==3) cout<<"proba after scaling="<<proba<<endl;
	}
      }
      else if(killEdges && residual >= upperLimit){
	//continue;
	proba=stimePDFQueue[pmtType]->Eval(upperLimit);
	if(pmtType==1) proba/=1;
	if(verbose==3 && pmtType==1) cout<<"PMT type = "<< pmtType <<", Upper limit = "<<upperLimit<<", residual = "<<residual<<", proba = "<<proba<<endl;
      }
      else if(killEdges && residual <= lowerLimit){
	//continue;
	proba=stimePDFQueue[pmtType]->Eval(lowerLimit);
	if(pmtType==1) proba/=1;
	if(verbose==3 && pmtType==1) cout<<"PMT type = "<< pmtType <<", Lower limit = "<<lowerLimit<<", residual = "<<residual<<", proba = "<<proba<<endl;
      }
      else{
	proba=0;
	//cout<<"GUAAAAAAAAAAAAAAAAAAAAAAAH, kill edges="<<killEdges<<", residual="<<residual<<endl;
      }
    }
    else{
      bool condition;
      if(pmtType == 0 || (pmtType == 1 && limitmPMT) ) condition = residual > hitTimeLimitsNegative && residual < hitTimeLimitsPositive;
      else condition = true;
      if(condition){
	NLL++;//= stimePDFQueue[pmtType]->Eval(residual);
	if(verbose==3){
	  cout<<"Residual="<<residual<<", pmt type="<<pmtType<<endl;
	  cout << "hit#"<<ihit<<", hit time =" << hitTime << ", vertex time = " << vertexPosition[0] << ", distance PMT vs vertex = " << distance << ", tof="<<tof<<endl;
	  cout<<"residual="<<residual<<", proba="<<proba<<endl;
	}
      }
      else if(killEdges) continue;
      else proba=0;
    }

    if(likelihood){
      if(proba<0){
	cout<<"Error in PDF"<<endl;
	if(proba>-1e-1) proba=0;//Since spline sometimes slightly goes below 0 due to interpolation
	else{cout <<"Error in "<< residual << "ns where proba = " << proba << endl;
	  return 0;
	}
      }
      if(proba==0) proba=1e-20;
      NLL += -TMath::Log(proba);
    }
  }
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
      NLLdir2=findNLLDirectionality(vertexPosition, nhits, verbose,lowerLimit,upperLimit);
    }
    else{
      //NLLdir=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,hitTimeLimitsNegative,hitTimeLimitsPositive);
      NLLdir2=findNLLDirectionality(vertexPosition, nhits, verbose,hitTimeLimitsNegative,hitTimeLimitsPositive);
    }
    if(verbose == 2) cout<<"NLL = "<<NLL<<", dir (L) = "<<TMath::Exp(-NLLdir)<<", dir 2 = "<<NLLdir2<<endl;
    if(directionality == 1) NLL+=NLLdir2;
    else if(directionality == 2) NLL=NLLdir2;
  }

  return NLL;
}

//Coarse search of the vertex. It will search in a cylinder having the radius = tankRadius and a height tankHeight.
//The step between each grid points is given by stepSize (in centimeters).
//We give the list of hit information through the 2D array essentialHitInfo to make the code faster, instead of using the whole code information
//The code provide as an output not only one vertex, but a list of vertices whose likelihood is higher (Negative Log Likelihood is lower). The number of output vertices can be chosen using "tolerance".
double ** searchVertex(double tankRadius,double tankHeight,double integrationTimeWindow,double stepSize,int nhits,int pmtType,double * trueVertexPosition,int tolerance=1,int verbose=0,bool likelihood=false,double lowerLimit=stimePDFLimitsQueueNegative, double upperLimit=stimePDFLimitsQueuePositive,int directionality=false){
  //2. How to set the search?
  //double stepSize = 0.5;//in m
  double * bestReconstructedVertexPosition = new double[4];
  
  vector < double[4] > reconstructedVertexPositionArray;
  reconstructedVertexPositionArray.clear();
  //double ** reconstructedVertexPosition = new double[4];
  double minNLL=1e15;
  std::map< double,double* > vertexContainer;
  //std::map< double,std::vector<double> > vertexContainer;
  clock_t timeStart=clock();
  
  //for(double time = -50; time < integrationTimeWindow; time+=(stepSize/lightSpeed)){
  for(double time = -50; time < 50; time+=(stepSize/lightSpeed)){
    if(verbose) cout << "Time = " << time << "ns" << endl;
    for(double radius = 0;radius <= tankRadius;radius+=stepSize){
      if(verbose) cout << "Radius = " << radius << "m" << endl;
      double perimeter = 2*TMath::Pi()*radius;
      int numberOfStepsOnCircle = floor(perimeter/stepSize);
      if(numberOfStepsOnCircle == 0) numberOfStepsOnCircle = 1;
      double angleStepOnCircle = 2*TMath::Pi()/numberOfStepsOnCircle;
      
      for(double angle=0;angle<=2*TMath::Pi();angle+=angleStepOnCircle){
	//if(verbose) cout << "Angle = " << angle << "rad" << endl;
	
	for(double height=-tankHeight/2.;height <= tankHeight/2.;height+=stepSize){
	  //if(verbose) cout << "Height = " << height << "m" << endl;
	  
	  double vertexPosition[4];//In centimeters
	  vertexPosition[0] = time;
	  vertexPosition[1] = radius*TMath::Cos(angle);
	  vertexPosition[2] = radius*TMath::Sin(angle);
	  vertexPosition[3] = height;

	  double distanceTrue=0;
	  for(int i=0;i<3;i++){
	    distanceTrue += pow(trueVertexPosition[i+1]-vertexPosition[i+1],2);
	  }
	  distanceTrue = TMath::Sqrt(distanceTrue);
	  if(verbose == 2) cout<<"time distance to true = "<<vertexPosition[0]-trueVertexPosition[0]<<"ns, distance to true = "<<distanceTrue<<"m, "<<endl;
	  //double NLL=findNLL(vertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit);
	  double NLL=findNLL(vertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);

	  NLL+=trand->Uniform(0,1e-5);//Here only not to have exactly the same value for 2NLL, and be removed from map
	  
	  //std::vector<double> candidateReconstructedVertexPosition;//[4];// = new double[4];
	  double * candidateReconstructedVertexPosition = new double[4];
	  for(int i=0;i<4;i++) candidateReconstructedVertexPosition[i] = vertexPosition[i];
	  //for(int i=0;i<4;i++) candidateReconstructedVertexPosition.push_back(vertexPosition[i]);
	  //reconstructedVertexPositionArray.push_back(candidateReconstructedVertexPosition)//;
	  //std::map< double, vector<double> >::iterator itMap = vertexContainer.end()-1;

	  vertexContainer[NLL] = candidateReconstructedVertexPosition;
	  if(vertexContainer.size()>tolerance){
	    vertexContainer.erase(--vertexContainer.rbegin().base());
	  }
	 
	  
	  if(NLL<minNLL){
	    minNLL=NLL;
	    //cout << "minNLL = "<<minNLL << endl;
	    for(int i=0;i<4;i++) bestReconstructedVertexPosition[i] = vertexPosition[i];
	  }
	  //candidateReconstructedVertexPosition.clear();
	  //if(verbose && distanceToTrue < 2) cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<endl;
	}
      }
    }
  }
  clock_t timeEndLoop=clock();
  double ** reconstructedVertexPosition = new double*[tolerance];
  for(int count=0;count<tolerance;count++){
    reconstructedVertexPosition[count] = new double[5];
    for(int i=0;i<5;i++) reconstructedVertexPosition[count][i]=0;
  }
  int counter=0;
  for(std::map< double, double * >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
  //for(std::map< double, vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
    if(counter<tolerance){
      //reconstructedVertexPosition[counter] = new double[4];
      for(int i=0;i<4;i++) reconstructedVertexPosition[counter][i] = (it->second)[i];
      reconstructedVertexPosition[counter][4] = it->first;
      //if(verbose) cout<<"NLL = "<<it->first<<", vertex = "<<it->second[0]<<", "<<it->second[1]<<", "<<it->second[2]<<", "<<it->second[3]<<endl;

      if(verbose){
	double distanceToTrue=0;
	for(int i=0;i<3;i++){
	  distanceToTrue += pow(trueVertexPosition[i+1]-it->second[i+1],2);
	}
	distanceToTrue = TMath::Sqrt(distanceToTrue);
	cout<<"NLL = "<<it->first<<", time distance to true = "<<it->second[0]-trueVertexPosition[0]<<"ns, distance to true = "<<distanceToTrue<<"m, "<<endl;
	if(verbose) cout<<"NLL = "<<it->first<<", vertex = "<<it->second[0]<<", "<<it->second[1]<<", "<<it->second[2]<<", "<<it->second[3]<<endl;
      }
    }
    counter++;
  }
  double distanceToTrue=0;
  for(int i=0;i<3;i++){
    distanceToTrue += pow(trueVertexPosition[i+1]-bestReconstructedVertexPosition[i+1],2);
  }
  distanceToTrue = TMath::Sqrt(distanceToTrue);
  if(verbose){
    cout << "Candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << bestReconstructedVertexPosition[0] - trueVertexPosition[0] << "ns, with a min NLL = " <<minNLL << endl;
    cout << "Time of the true event = " << trueVertexPosition[0] << endl;
    }
  clock_t timeEnd=clock();
  cout<<"Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<endl<<endl;
  return reconstructedVertexPosition;
  //return bestReconstructedVertexPosition;

  for(int i=0;i<tolerance;i++) delete reconstructedVertexPosition[i];
  delete reconstructedVertexPosition;
  delete bestReconstructedVertexPosition;
  //delete candidateReconstructedVertexPosition;
}


//Fine grid search of the vertex. It will search in a box this time, and not a cylinder as searchVertex is doing. It is useful when we have one or several first candidate vertices. Therefore, this code is generally run after serachVertex.
//The box is centered around the initial vertex provided.
//The box size is 2 x limits in each direction, centered on the initial vertex provided.
//nCandidate provide the size of the list of inital vertices, while tolerance provide the size of the output list.
//Other informations are the same as searchVertex.
double ** searchVertexFine(double ** initialVertex,double * limits, double stepSize,int nhits,int pmtType,double * trueVertexPosition,int nCandidates = 1,int tolerance = 1,int verbose=0,bool likelihood=false,bool average=false,double lowerLimit=stimePDFLimitsQueueNegative, double upperLimit=stimePDFLimitsQueuePositive, int directionality=false){
  cout<<"Start fine search, use directionality? "<<directionality<<endl;
  //2. How to set the search?
  //double stepSize = 0.5;//in m
  //double * reconstructedVertexPosition = new double[4];
  double * bestReconstructedVertexPosition = new double[4];
  double minNLL=1e15;
  clock_t timeStart=clock();

  //std::map< double,std::vector<double> > vertexContainer;
  std::map< double,double* > vertexContainer;
  
  for(int icand=0;icand<nCandidates;icand++){
    cout<<"Candidate #"<<icand<<endl;
    for(double time = initialVertex[icand][0] - limits[0]/lightSpeed; time <= initialVertex[icand][0]+limits[0]/lightSpeed; time+=(stepSize/lightSpeed)){
      if(verbose == 2) cout << "Time = " << time << "ns" << endl;
      for(double x = initialVertex[icand][1] - limits[1]; x <= initialVertex[icand][1] + limits[1];x+=stepSize){
	if(verbose == 2) cout << "x = " << x << "m" << endl;
	for(double y = initialVertex[icand][2] - limits[2]; y <= initialVertex[icand][2] + limits[2];y+=stepSize){
	  if(verbose == 2) cout << "y = " << y << "m" << endl;
	  for(double z = initialVertex[icand][3] - limits[3]; z <= initialVertex[icand][3] + limits[3];z+=stepSize){
	    if(verbose == 2) cout << "z = " << z << "m" << endl;

	    double vertexPosition[4];//In centimeters
	    vertexPosition[0] = time;
	    vertexPosition[1] = x;//radius*TMath::Cos(angle);
	    vertexPosition[2] = y;//radius*TMath::Sin(angle);
	    vertexPosition[3] = z;//height;

	    double distanceTrue=0;
	    for(int i=0;i<3;i++){
	      distanceTrue += pow(trueVertexPosition[i+1]-vertexPosition[i+1],2);
	    }
	    distanceTrue = TMath::Sqrt(distanceTrue);
	    if(verbose == 2) cout<<"time distance to true = "<<vertexPosition[0]-trueVertexPosition[0]<<"ns, distance to true = "<<distanceTrue<<"m, "<<endl;

	    double NLL=0;
	    if(average){
	      if(verbose==2) cout << endl << "Vertex timing = "<<vertexPosition[0]-trueVertexPosition[0]<<", (x,y,z)=("<<vertexPosition[1]-trueVertexPosition[1]<<","<<vertexPosition[2]-trueVertexPosition[2]<<","<<vertexPosition[3]-trueVertexPosition[3]<<")"<<endl;
	      //double vNLL[nAveraging];
	      for(int irand=0;irand<nAveraging;irand++){
		//Determine the position around the vertex
		double radius = 1.5*trand->Uniform(0,stepSize);
		double phi = trand->Uniform(0,2*TMath::Pi());
		double theta = trand->Uniform(0,TMath::Pi());
		double timeRand = 1.5*trand->Uniform(-stepSize/lightSpeed,stepSize/lightSpeed);
		double randVertexPosition[4];//In centimeters
		randVertexPosition[0] = vertexPosition[0] + timeRand;//rand->Uniform(-stepSize/lightSpeed/2,stepSize/lightSpeed/2);
		randVertexPosition[1] = vertexPosition[1] + radius*TMath::Sin(theta)*TMath::Cos(phi);
		randVertexPosition[2] = vertexPosition[2] + radius*TMath::Sin(theta)*TMath::Sin(phi);
		randVertexPosition[3] = vertexPosition[3] + radius*TMath::Cos(theta);
		//Evaluate its NLL
		double nll=findNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);
		//double nll=findNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit);
		//double nll=findRMS(pmtType,randVertexPosition,essentialHitInfo,nhits,verbose,lowerLimit,upperLimit);
		if(verbose==2) cout << "Randomized Vertex timing = "<<randVertexPosition[0]-trueVertexPosition[0]<<", (x,y,z)=("<<randVertexPosition[1]-trueVertexPosition[1]<<","<<randVertexPosition[2]-trueVertexPosition[2]<<","<<randVertexPosition[3]-trueVertexPosition[3]<<"), nll="<<nll<<endl;
		NLL+=nll;
		//vNLL[irand]=nll;
	      }
	      
	      if(nAveraging!=0) NLL/=nAveraging;
	      //double rmsNLL=0;
	      //for(int irand=0;irand<nAveraging;irand++){
	      //rmsNLL+=pow(vNLL[irand]-NLL,2);
	      //}
	      //if(nAveraging!=0) rmsNLL/=nAveraging;
	      //rmsNLL=TMath::Sqrt(rmsNLL);
	      //NLL=rmsNLL;
	    }
	    else NLL=findNLL(vertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit,false,false,directionality);
	    //else NLL=findChi2BinnedExpo(pmtType,vertexPosition,nhits,1,100);
	    //else NLL=findRMS(pmtType,vertexPosition,essentialHitInfo,nhits,verbose,lowerLimit,upperLimit);

	    
	    double * candidateReconstructedVertexPosition = new double[4];
	    for(int i=0;i<4;i++) candidateReconstructedVertexPosition[i] = vertexPosition[i];
	    //for(int i=0;i<4;i++) candidateReconstructedVertexPosition.push_back(vertexPosition[i]);
	    //reconstructedVertexPositionArray.push_back(candidateReconstructedVertexPosition)//;
	    //std::map< double, vector<double> >::iterator itMap = vertexContainer.end()-1;
	    
	    vertexContainer[NLL] = candidateReconstructedVertexPosition;
	    if(vertexContainer.size()>tolerance){
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
	    //cout << "minNLL = "<<minNLL << endl;
	    for(int i=0;i<4;i++) bestReconstructedVertexPosition[i] = vertexPosition[i];
	  }/*
	  candidateReconstructedVertexPosition.clear();*/
	  //if(verbose && distanceToTrue < 2) cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<endl;
	  //if(verbose && y=) cout << "Distance to true vertex = " << distanceToTrue << ", Probability = " << NLL <<endl;
	}
      }
      }
    }
  }
  clock_t timeEndLoop=clock();
  double ** reconstructedVertexPosition = new double*[tolerance];
  for(int count=0;count<tolerance;count++){
    reconstructedVertexPosition[count] = new double[5];
    for(int i=0;i<5;i++) reconstructedVertexPosition[count][i]=0;
  }
  int counter=0;
  for(std::map< double, double* >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
    //for(std::map< double, vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
    if(counter<tolerance){
      //reconstructedVertexPosition[counter] = new double[4];
      for(int i=0;i<4;i++) reconstructedVertexPosition[counter][i] = (it->second)[i];
      reconstructedVertexPosition[counter][4] = it->first;

      if(verbose){
	double distanceToTrue=0;
	for(int i=0;i<3;i++){
	  distanceToTrue += pow(trueVertexPosition[i+1]-it->second[i+1],2);
	}
	distanceToTrue = TMath::Sqrt(distanceToTrue);
	cout<<"NLL = "<<it->first<<", time distance to true = "<<it->second[0]-trueVertexPosition[0]<<"ns, distance to true = "<<distanceToTrue<<"m, "<<endl;
      //if(verbose) cout<<"NLL = "<<it->first<<", vertex = "<<it->second[0]<<", "<<it->second[1]<<", "<<it->second[2]<<", "<<it->second[3]<<endl;
      }
    }
    
    counter++;
  }
  double distanceToTrue=0;
  for(int i=0;i<3;i++){
    distanceToTrue += pow(trueVertexPosition[i+1]-bestReconstructedVertexPosition[i+1],2);
  }
  distanceToTrue = TMath::Sqrt(distanceToTrue);
  if(verbose){
    cout << "Fine candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << bestReconstructedVertexPosition[0] - trueVertexPosition[0] << "ns, with a min NLL = " <<minNLL << endl;
    cout << "Time of the true event = " << trueVertexPosition[0] << endl;
  }
  clock_t timeEnd=clock();
  cout<<"Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<endl<<endl;
  
  return reconstructedVertexPosition;
  //return bestReconstructedVertexPosition;
  
  for(int i=0;i<tolerance;i++) delete reconstructedVertexPosition[i];
  delete reconstructedVertexPosition;
  delete bestReconstructedVertexPosition;
  //delete candidateReconstructedVertexPosition;
}



void minuitLikelihood(int& nDim, double * gout, double & NLL, double par[], int flg){

  double vertexPosition[4];//In centimeters
  for(int i=0;i<4;i++) vertexPosition[i]=par[i];
  int pmtType=par[4];
  int nhits=par[5];
  double lowerLimit=par[6];
  double upperLimit=par[7];
  double expoSigma=par[8];
  double directionality=par[9];
  //double scalingFactor=par[8];
  //NLL=findNLLSignalXDR(pmtType,vertexPosition,nhits,true,1,lowerLimit,upperLimit);
  //NLL=findNLLExpo(pmtType,vertexPosition,nhits,true,1,lowerLimit,upperLimit);
  //NLL=findNLL(vertexPosition,nhits,true,1,lowerLimit,upperLimit,true,false,directionality);
  NLL=findNLL(vertexPosition,nhits,true,1,lowerLimit,upperLimit,true,false,directionality);
  //cout<<"Directionality"<<endl;
  //NLL=findRateAndShape(pmtType,vertexPosition,nhits,1,lowerLimit,upperLimit,scalingFactor,true);
  //NLL=findChi2Binned(pmtType,vertexPosition,nhits,1,lowerLimit,upperLimit);
  //NLL=findChi2BinnedExpo(pmtType,vertexPosition,nhits,1,expoSigma);
  //if(verbose){
  //cout<<setprecision(10)<<"NLL="<<NLL<<", Vertex: "<<vertexPosition[0]<<"ns, ("<<vertexPosition[1]<<", "<<vertexPosition[2]<<", "<<vertexPosition[3]<<")"<<endl;
  //cout<<"Expo sigma = "<<expoSigma<<"ns"<<endl;
  //}
  //cout<<"Timing window = "<<lowerLimit<<"ns - "<<upperLimit<<"ns"<<endl<<endl;
  //double nll=findNLL(randVertexPosition,nhits,likelihood,verbose,lowerLimit,upperLimit);

}

//here is the function where minimization gonna take place. The loop is inside this function
double ** minimizeVertex(double ** initialVertex,double * limits, double stepSize,int nhits,int pmtType,double * trueVertexPosition,int nCandidates = 1,int tolerance = 1,int verbose=0,bool likelihood=false,bool average=false,double lowerLimit=stimePDFLimitsQueueNegative, double upperLimit=stimePDFLimitsQueuePositive, int directionality = true){
  cout<<"Minimizer"<<endl;
  //2. How to set the search?
  //double stepSize = 0.5;//in m
  //double * reconstructedVertexPosition = new double[4];
  double * bestReconstructedVertexPosition = new double[4];
  double minNLL=1e15;
  clock_t timeStart=clock();

  //std::map< double,std::vector<double> > vertexContainer;
  std::map< double,double* > vertexContainer;

  if(verbose==2) cout<<"Minimizer"<<endl;
  TFitter * minimizer = new TFitter(8);//4=nb de params?
  TMinuit * minuit = minimizer->GetMinuit();
  double arglist[20];
  int err=0;
  //arglist[0]=0;
  double p1=verbose-1;
  minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
  minuit->SetErrorDef(1);
  arglist[0]=2;
  minuit->mnexcm("SET STR",arglist,1,err);//set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
  minimizer->SetFCN(minuitLikelihood);//here is function to minimize

  double Mig[2]={1e6,1e0};//maxcalls and tolerance
  double Imp2[1]={10000};
  for(int i=0;i<4;i++){
    minimizer->SetParameter(i,Form("vertex%d",i),0,stepSize/(i==0?lightSpeed:1),-limits[i]/(i==0?lightSpeed:1),limits[i]/(i==0?lightSpeed:1));
  }
  minimizer->SetParameter(4,"pmtType",pmtType,pmtType,pmtType-1,pmtType+1);
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
    if(verbose>=1) cout<<"Candidate #"<<icand<<endl;
      for(int i=0;i<4;i++){
	//initialVertex[icand][i]-=trueVertexPosition[i];
	cout<<"Initial="<<initialVertex[icand][i]<<endl;
	minimizer->SetParameter(i,Form("vertex%d",i),initialVertex[icand][i],stepSize/(i==0?lightSpeed:1),initialVertex[icand][i]-limits[i]/(i==0?lightSpeed:1),initialVertex[icand][i]+limits[i]/(i==0?lightSpeed:1));
      }
      minimizer->SetParameter(4,"pmtType",pmtType,pmtType,pmtType-1,pmtType+1);
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
      
      //double nll=findChi2BinnedExpo(pmtType,initialVertex[icand],nhits,1,100);
      //double nll=findNLL(initialVertex[icand],nhits,likelihood,2,lowerLimit,upperLimit,true,true);
      //double nll2=findNLL(initialVertex[icand],nhits,likelihood,2,-3,4);
      //cout<<"NLL="<<nll<<endl;
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
      //double nll=findChi2Binned(pmtType,initialVertex[icand],nhits,1,lowerLimit,upperLimit);
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


      double * candidateReconstructedVertexPosition = new double[4];
      for(int i=0;i<4;i++) candidateReconstructedVertexPosition[i] = minimizer->GetParameter(i);
      minuit=minimizer->GetMinuit();
      double fmin, fedm, errdef;
      int npar,nparx,istat;
      minuit->mnstat(fmin,fedm,errdef,npar,nparx,istat);
      double NLL=fmin;
      if(verbose>=1) cout<<"NLL="<<NLL<<endl;
      vertexContainer[NLL] = candidateReconstructedVertexPosition;
      if(vertexContainer.size()>tolerance){
	vertexContainer.erase(--vertexContainer.rbegin().base());
      }

      if(NLL<minNLL){
	minNLL=NLL;
	//cout << "minNLL = "<<minNLL << endl;
	for(int i=0;i<4;i++) bestReconstructedVertexPosition[i] = candidateReconstructedVertexPosition[i];
      }
  }	 
	 
  clock_t timeEndLoop=clock();
  double ** reconstructedVertexPosition = new double*[tolerance];
  for(int count=0;count<tolerance;count++){
    reconstructedVertexPosition[count] = new double[5];
    for(int i=0;i<5;i++) reconstructedVertexPosition[count][i]=0;
  }
  int counter=0;
  for(std::map< double, double* >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
    //for(std::map< double, vector<double> >::iterator it = vertexContainer.begin();it!=vertexContainer.end();it++){
	if(counter<tolerance){
	  //reconstructedVertexPosition[counter] = new double[4];
	  for(int i=0;i<4;i++) reconstructedVertexPosition[counter][i] = (it->second)[i];
	  reconstructedVertexPosition[counter][4] = it->first;
	  
	  if(verbose){
	    double distanceToTrue=0;
	    for(int i=0;i<3;i++){
	      distanceToTrue += pow(trueVertexPosition[i+1]-it->second[i+1],2);
	    }
	    distanceToTrue = TMath::Sqrt(distanceToTrue);
	    cout<<"NLL = "<<it->first<<", time distance to true = "<<it->second[0]-trueVertexPosition[0]<<"ns, distance to true = "<<distanceToTrue<<"m, "<<endl;
	    //if(verbose) cout<<"NLL = "<<it->first<<", vertex = "<<it->second[0]<<", "<<it->second[1]<<", "<<it->second[2]<<", "<<it->second[3]<<endl;
	  }
	}
	
	counter++;
  }
  double distanceToTrue=0;
  for(int i=0;i<3;i++){
    distanceToTrue += pow(trueVertexPosition[i+1]-bestReconstructedVertexPosition[i+1],2);
  }
  distanceToTrue = TMath::Sqrt(distanceToTrue);
  if(verbose){
    cout << "Fine candidate is at a distance " << distanceToTrue << "m" << ", time difference = " << bestReconstructedVertexPosition[0] - trueVertexPosition[0] << "ns, with a min NLL = " <<minNLL << endl;
    cout << "Time of the true event = " << trueVertexPosition[0] << endl;
  }
  clock_t timeEnd=clock();
  cout<<"Time it took for the loop = "<<(timeEndLoop - timeStart)/1e6<<", time total = "<<(timeEnd - timeStart)/1e6<<endl<<endl;

  return reconstructedVertexPosition;
  //return bestReconstructedVertexPosition;
  
  for(int i=0;i<tolerance;i++) delete reconstructedVertexPosition[i];
  delete reconstructedVertexPosition;
  delete bestReconstructedVertexPosition;
  //delete candidateReconstructedVertexPosition;

}

int main (int argc, char **argv){

  loadSplines();
  
  std::string filename="";
  std::string outfilename="";
  int verbose=-1;
  bool hybrid = true;//false;//true;
  bool Gamma = false;//Is the mother particle a gamma or another particle?
  bool plotDigitized = true;//false;//true;//false;//true;//false;//true;//false;//true;//false;

  bool HK=true;
  double TotalNHits = 1e4;
  double TankSize = 1100;//in cm, the maximal size of the plots for vertices etc...
  double TankRadius = 742/2;//in cm, the maximal size of the plots for vertices etc...
  double TankHalfHeight = 1042/2;//in cm, the maximal size of the plots for vertices etc...

  int startEvent=0;
  int endEvent=0;
  int c=-1;
  while( (c = getopt(argc,argv,"f:o:s:e:h")) != -1 ){//input in c the argument (-f etc...) and in optarg the next argument. When the above test becomes -1, it means it fails to find a new argument.
      switch(c){
      case 'f':
	filename = optarg;
	break;
      case 'h':
	HK = true;
	break;
      case 'o':
	outfilename = optarg;
	break;
      case 's':
	startEvent = atoi(optarg);
	break;
      case 'e':
	endEvent = atoi(optarg);
	break;
      }
  }
  
  if(HK){
    TotalNHits = 3e4;
    TankSize = 7100;//in cm, the maximal size of the plots for vertices etc...
    TankRadius = 7080/2;//in cm, the maximal size of the plots for vertices etc...
    TankHalfHeight = 5480/2;//in cm, the maximal size of the plots for vertices etc...
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );
  }
    //nevt = std::min( 1e3 , ((double) WSTreeHits->GetEntries()) );

  
  TFile *file;
  // Open the file
  if (filename == ""){
    // Test input
    file = new TFile("/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_e10_center_nominal_fulltank_0hitstrigger_10000.root","read");


    
  }else{
    file = new TFile(filename.c_str(),"read");
  }
  if (!file->IsOpen()){
    cout << "Error, could not open input file: " << filename << endl;
    return -1;
  }
  
  // Get the a pointer to the tree from the file
  TTree *tree = (TTree*)file->Get("wcsimT");
  
  // Get the number of events
  int nevent = std::min(((int)tree->GetEntries()),20000);
  if(endEvent!=0) nevent = endEvent;
  if(verbose) printf("nevent %d\n",nevent);
  
  // Create a WCSimRootEvent to put stuff from the tree in

  WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
  WCSimRootEvent* wcsimrootsuperevent2 = new WCSimRootEvent();

  // Set the branch address for reading from the tree
  TBranch *branch = tree->GetBranch("wcsimrootevent");
  branch->SetAddress(&wcsimrootsuperevent);
  // Force deletion to prevent memory leak 
  tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

  TBranch *branch2;
  if(hybrid){
    branch2 = tree->GetBranch("wcsimrootevent2");
    branch2->SetAddress(&wcsimrootsuperevent2);
  // Force deletion to prevent memory leak 
    tree->GetBranch("wcsimrootevent2")->SetAutoDelete(kTRUE);
  }

  // Geometry tree - only need 1 "event"
  TTree *geotree = (TTree*)file->Get("wcsimGeoT");
  geotree->SetBranchAddress("wcsimrootgeom", &geo);
  if(verbose) std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
  if (geotree->GetEntries() == 0) {
      exit(9);
  }
  geotree->GetEntry(0);
  cout << "Number of PMTs of 1st type = " << geo->GetWCNumPMT() << endl;
  cout << "Number of PMTs of 2nd type = " << geo->GetWCNumPMT(true) << endl;


  int nMasked=0;
  bool maskedTube[nPMTtypes][maxPMTTubes]={{false}};
  if(maskPMT){
    trand2->SetSeed(1);
    int numPMT=(pmtTypeMasked==0?geo->GetWCNumPMT():geo->GetWCNumPMT(true));
    for (int icab=0; icab<numPMT; icab++) { 
     int tubeNumber=icab+1;
     //pmt = geo->GetPMT(tubeNumber-1,pmtTypeMasked==0?false:true);//check icab -1 or icab
      double maskingRand=trand2->Uniform();
      if(maskingRand<=maskedProbability){
	//cout<<"masked is icab = "<<icab<<endl;
	maskedTube[pmtTypeMasked][tubeNumber]=true;
	nMasked++;
      }
    }
  }
  cout<<"number of PMTs masked = "<<nMasked<<endl;
  
  // Options tree - only need 1 "event"
  TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
  WCSimRootOptions *opt = 0; 
  opttree->SetBranchAddress("wcsimrootoptions", &opt);
  if(verbose) std::cout << "Optree has " << opttree->GetEntries() << " entries" << std::endl;
  if (opttree->GetEntries() == 0) {
    exit(9);
  }
  opttree->GetEntry(0);
  opt->Print();

  // start with the main "subevent", as it contains most of the info
  // and always exists.
  WCSimRootTrigger* wcsimrootevent;
  WCSimRootTrigger* wcsimrootevent2;


  Float_t bsvertex[4],bsresult[6], bsdir[3]; //*cm
  Float_t bsgood[10];
  Int_t bsnhit[2]; //nsel (SLE)
  Int_t bsnhit_pmtType0[2]; //Only for firs PMT type. Allow us to do a fit using only first PMT type hits
  Int_t bsnsel[2]; //nsel (SLE)
  float mcvertex[4], mcdir[3], mcE, mcM; //*cm ...?
  float mcTime;
  float *chi2[nPMTtypes];
  float goodness[nPMTtypes];
  double timeWindow[nPMTtypes][2];
  Int_t bsnhit_intime[nPMTtypes]; //Hit in the time window considered
  Int_t bsnhit_200ns[nPMTtypes]; //nsel (SLE)
  Int_t pmtTypeFill;
  
  if(outfilename == "") outfilename = "out.root";
  TFile * outfile = new TFile(outfilename.c_str(),"RECREATE");
  TTree *bstree = new TTree("bstree","bstree");
  bstree->Branch("bsvertex", bsvertex, "bsvertex[4]/F");
  bstree->Branch("mcvertex", mcvertex, "mcvertex[4]/F");
  bstree->Branch("mcM", &mcM, "mcM/F");
  bstree->Branch("mcE", &mcE, "mcE/F");
  bstree->Branch("bsdir", bsdir, "bsdir[3]/F");
  bstree->Branch("mcdir", mcdir, "mcdir[3]/F");
  bstree->Branch("bsgoodness", &bsgood[2], "bsgoodness/F");
  bstree->Branch("bsnhit", bsnhit, "bsnhit[2]/I");
  bstree->Branch("bsnhit_pmtType0", bsnhit_pmtType0, "bsnhit_pmtType0[2]/I");
  bstree->Branch("bsgoodness_BQ", goodness, Form("bsgoodness_BQ[%d]/F",nPMTtypes));
  bstree->Branch("bsnhit_intime", bsnhit_intime, Form("bsnhit_intime[%d]/I",nPMTtypes));
  bstree->Branch("bsnhit_200ns", bsnhit_200ns, Form("bsnhit_200ns[%d]/I",nPMTtypes));
  bstree->Branch("pmtType", &pmtTypeFill, "pmtType/I");  //B.Q 
  
  TH1F *h1 = new TH1F("PMT Hits", "PMT Hits", 8000, 0, 8000);
  TH1F *hvtx0 = new TH1F("Event VTX0", "Event VTX0", 200, -1500, 1500);
  TH1F *hvtx1 = new TH1F("Event VTX1", "Event VTX1", 200, -1500, 1500);
  TH1F *hvtx2 = new TH1F("Event VTX2", "Event VTX2", 200, -1500, 1500);
  
  int num_trig=0;

    //TFile * fOutput = new TFile(OutputFile,"recreate");

  for(int i=0;i<nPMTtypes;i++){
    std::cout << " PMT Type: " << i << std::endl;
    ChargeProfile2D_onlyFront[i] = new TH2D(Form("ChargeProfile2D_onlyFront_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2D[i] = new TH2D(Form("ChargeProfile2D_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2DTop[i] = new TH2D(Form("ChargeProfile2DTop_pmtType%d",i),"",TankSize/2,-TankSize/2.,TankSize/2.,TankSize/2,-TankSize/2.,TankSize/2.);
    //ChargeProfile2DCylindricCoordinate[i] = new TH2D(Form("ChargeProfile2DCylindricCoordinate_pmtType%d",i),"",3600,0,360,TankSize,-TankSize,TankSize);
    ChargeProfile2DCylindricCoordinate[i] = new TH2D(Form("ChargeProfile2DCylindricCoordinate_pmtType%d",i),"",720,0,360,TankSize/2,-TankSize/2.,TankSize/2.);
    ChargeProfile2DRelativeAngle[i] = new TH2D(Form("ChargeProfile2DRelativeAngle_pmtType%d",i),"",720,0,360,TankSize/2,-TankSize/2.,TankSize/2.);
    TimeTOFProfileXTOF[i] = new TH2D(Form("TimeTOFProfileXTOF_pmtType%d",i),"",1e4,-1e2,1e3,1e2,0,5e2);TimeTOFProfileXTOF[i]->Sumw2();
    
    ChargeProfile[i] = new TH1D(Form("ChargeProfile_pmtType%d",i),"",720,0,180);
    HitProfile[i] = new TH1D(Form("HitProfile_pmtType%d",i),"",720,0,180);

    TimeProfile[i] = new TH1D(Form("TimeProfile_pmtType%d",i),"",1e4,0,5e3);
    TimeHitProfile[i] = new TH1D(Form("TimeHitProfile_pmtType%d",i),"",1e4,0,5e3);

    TimeTOFProfile[i] = new TH1D(Form("TimeTOFProfile_pmtType%d",i),"",1e4,-1e2,1e3);
    TimeTOFProfile[i]->Sumw2();
    TimeTOFProfile_currentEvent[i] = new TH1D(Form("TimeTOFProfile_currentEvent_pmtType%d",i),"",1e4,-1e2,1e3);
    TimeTOFProfile[i]->Sumw2();
    HitTimeTOFProfile[i] = new TH1D(Form("HitTimeTOFProfile_pmtType%d",i),"",1e4,-1e2,1e3);
    HitTimeTOFProfile[i]->Sumw2();

    ChargePerPMT[i] = new TH1D(Form("ChargePerPMT_pmtType%d",i),"",500,0,500);

    TotalCharge[i] = new TH1D(Form("TotalCharge_pmtType%d",i),"",1e4,0,3e4);
    TotalHit[i] = new TH1D(Form("TotalHit_pmtType%d",i),"",1e4,0,TotalNHits);

    TotalChargeXdWall[i] = new TH2D(Form("TotalChargeXdWall_pmtType%d",i),"",TankSize/2,0,TankSize/2.,1e4,0,3e4);

    //EfficiencyLEXdWall[i] = new TH1D(Form("EfficiencyLEXdWall_pmtType%d",i),"",TankSize/2,0,TankSize/2.,1e4,0,3e4);
    //EfficiencyLEXdWallNorm[i] = new TH1D(Form("EfficiencyLEXdWallNorm_pmtType%d",i),"",TankSize/2,0,TankSize/2.,1e4,0,3e4);
    timeWindowXsignal[i] = new TH3D(Form("timeWindowXsignal_pmtType%d",i),"",30,0,300,20,-10,0,40,0,20);//p.e vs lower limit on time window vs upper limit.
    timeWindowXDR[i] = new TH3D(Form("timeWindowXDR_pmtType%d",i),"",30,0,300,20,-10,0,40,0,20);//p.e vs lower limit on time window vs upper limit.
    
    fitTimeToTrue[i] = new TH1D(Form("fitTimeToTrue_pmtType%d",i),"",1e4,-1000,1000);//in ns.
    fitDistanceToTrue[i] = new TH1D(Form("fitDistanceToTrue_pmtType%d",i),"",5e3,0,5e3);//in cm.
    fit4DdistanceXNLLdifference[i] = new TH2D(Form("fit4DdistanceXNLLdifference_pmtType%d",i),"",5e3,0,5e3,100,-1,1);

    fitNLLprofile_currentEvent[i] = new TH2D(Form("fitNLLprofile_currentEvent_pmtType%d",i),"",200,-10,10,400,0,200);
    normfitNLLprofile_currentEvent[i] = new TH2D(Form("normfitNLLprofile_currentEvent_pmtType%d",i),"",200,-10,10,400,0,200);
    fitNLLprofile_currentEvent[i]->Sumw2(); normfitNLLprofile_currentEvent[i]->Sumw2();

   
    ChargeProfile2D[i]->Sumw2();
    ChargeProfile2D_onlyFront[i]->Sumw2();
    ChargeProfile2DTop[i]->Sumw2();
    ChargeProfile2DCylindricCoordinate[i]->Sumw2();
    ChargeProfile2DRelativeAngle[i]->Sumw2();
    ChargeProfile[i]->Sumw2(); HitProfile[i]->Sumw2();TimeProfile[i]->Sumw2();TimeHitProfile[i]->Sumw2();ChargePerPMT[i]->Sumw2();
    TotalCharge[i]->Sumw2();TotalHit[i]->Sumw2();
    TotalChargeXdWall[i]->Sumw2();
  }
  /*
  
  TH2D * HitProfile2D = new TH2D("HitProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * TimeProfile2D = new TH2D("TimeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);
  TH2D * ChargeXPositionXTime = new TH2D("ChargeXPositionXTime","",36,0,180,TankSize,0,TankSize);
  TH2D * ChargeProfileXdWall = new TH2D("ChargeProfileXdWall","",720,0,180,100,0,TankSize);
  TH2D * TotalChargeXdWall = new TH2D("TotalChargeXdWall","",100,0,TankSize,1e3,0,1e5);
  TH2D * TimeAngleProfile = new TH2D("TimeAngleHitProfile","",720,0,180,1e4,0,5e3);TimeAngleProfile->Sumw2();
  TH2D * HitTimeAngleProfile = new TH2D("HitTimeAngleHitProfile","",720,0,180,1e4,0,5e3);HitTimeAngleProfile->Sumw2();
  TH2D * TimeTOFAngleProfile = new TH2D("TimeTOFAngleProfile","",720,0,180,1e4,-1e2,5e3);TimeTOFAngleProfile->Sumw2();
  TH2D * HitTimeTOFAngleProfile = new TH2D("HitTimeTOFAngleProfile","",720,0,180,1e4,-1e2,5e3);HitTimeTOFAngleProfile->Sumw2();
  TH2D * TimeTOFTOF = new TH2D("TimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);TimeTOFTOF->Sumw2();
  TH2D * HitTimeTOFTOF = new TH2D("HitTimeTOFTOF","",1e4,-1e2,5e3,1e2,0,5e2);HitTimeTOFTOF->Sumw2();
  TH3D * VertexPosition = new TH3D("VertexPosition","",100,-TankSize,TankSize,100,-TankSize,TankSize,100,-TankSize,TankSize);
  TH3D * VertexDirection = new TH3D("VertexDirection","",100,-1.1,1.1,100,-1.1,1.1,100,-1.1,1.1);
  TH1D * VertexXdWall = new TH1D("VertexXdWall","",100,0,TankSize);
  TH1D * ParentFlyingDistance = new TH1D("ParentFlyingDistance","",10000,0,TankSize);
  TH2D * EGammaSeparationXdWall = new TH2D("EGammaSeparationxdWall","",100,0,TankSize,1000,0.,5.);
  TH1D * EGammaSeparation = new TH1D("EGammaSeparation","",1000,0.,5.);
  //Set the 2D histogram to 0 to increase lisibility of plots
  for(int ibinx=1;ibinx<=ChargeProfile2D->GetNbinsX();ibinx++){
    for(int ibiny=1;ibiny<=ChargeProfile2D->GetNbinsY();ibiny++){
      HitProfile2D->SetBinContent(ibinx,ibiny,0);
      ChargeProfile2D->SetBinContent(ibinx,ibiny,0);     
      ChargeProfile2DXZ->SetBinContent(ibinx,ibiny,0);     
    }
  }
  //
  double ConversionLength = 36;//cm
  TH2D * ChargeProfile2D_PerEvent[nevt];
  */
  
  //TH2D * ChargeProfile2D = new TH2D("ChargeProfile2D","",TankSize/2,-TankSize,TankSize,TankSize/2,-TankSize,TankSize);

  // Now loop over events
  for (int ev=startEvent; ev<nevent; ev++)
  {
    // Read the event from the tree into the WCSimRootEvent instance
    tree->GetEntry(ev);      
    wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);
    if(hybrid) wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    //wcsimrootevent2 = wcsimrootsuperevent2->GetTrigger(0);
    if(verbose){
      printf("********************************************************");
      printf("Evt, date %d %d\n", wcsimrootevent->GetHeader()->GetEvtNum(),
	     wcsimrootevent->GetHeader()->GetDate());
      printf("Mode %d\n", wcsimrootevent->GetMode());
      printf("Number of subevents %d\n",
	     wcsimrootsuperevent->GetNumberOfSubEvents());
      
      printf("Vtxvol %d\n", wcsimrootevent->GetVtxvol());
      printf("Vtx %f %f %f\n", wcsimrootevent->GetVtx(0),
	     wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
    }
    hvtx0->Fill(wcsimrootevent->GetVtx(0));
    hvtx1->Fill(wcsimrootevent->GetVtx(1));
    hvtx2->Fill(wcsimrootevent->GetVtx(2));

    if(verbose){
      printf("Jmu %d\n", wcsimrootevent->GetJmu());
      printf("Npar %d\n", wcsimrootevent->GetNpar());
      printf("Ntrack %d\n", wcsimrootevent->GetNtrack());
      
    }

    std::vector<float> triggerInfo;
    triggerInfo.clear();
    triggerInfo = wcsimrootevent->GetTriggerInfo();

    std::vector<float> triggerInfo2;
    triggerInfo2.clear();
    if(hybrid) triggerInfo2 = wcsimrootevent2->GetTriggerInfo();

    if(verbose){
      for(int v=0;v<triggerInfo.size();v++){
	cout << "Trigger entry #" << v << ", info = " << triggerInfo[v] << endl;
      }
      if(hybrid){
	for(int v=0;v<triggerInfo2.size();v++){
	  cout << "Trigger2 entry #" << v << ", info = " << triggerInfo2[v] << endl;
	}
      }
    }

    double triggerShift[nPMTtypes];
    double triggerTime[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      triggerShift[pmtType]=0;
      triggerTime[pmtType]=0;
      if(triggerInfo.size()>=3){
	if(pmtType==0){
	  triggerShift[pmtType] = triggerInfo[1];
	  triggerTime[pmtType] = triggerInfo[2];
	}
      }
      if(triggerInfo2.size()>=3){
	if(pmtType==1 && hybrid){
	  triggerShift[pmtType] = triggerInfo2[1];
	  triggerTime[pmtType] = triggerInfo2[2];
	}
      }
    }    
    // Now read the tracks in the event
    
    // Get the number of tracks
    int ntrack = wcsimrootevent->GetNtrack();
    if(verbose) printf("ntracks=%d\n",ntrack);

    double particleStart[3];
    double particleStop[3];
    double particleDir[3];

    int i;
    // Loop through elements in the TClonesArray of WCSimTracks
    for (i=0; i<ntrack; i++)
    {
      TObject *element = (wcsimrootevent->GetTracks())->At(i);
      
      //WCSimRootTrack *wcsimroottrack = dynamic_cast<WCSimRootTrack*>(element);
      WCSimRootTrack *wcsimroottrack = (WCSimRootTrack*) (element);

      if(i==(ntrack-1)){//Mother particle
	for(int j=0; j<3; j++){
	  particleStart[j] = wcsimroottrack->GetStart(j);
	  particleStop[j] = wcsimroottrack->GetStop(j);
	  particleDir[j] = wcsimroottrack->GetDir(j);
	}
	//if(particleDir[0] != 1 || particleDir[1] != 0 || particleDir[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle direction has an issue" << endl << "###############################################" << endl;
	//if(particleStart[0] != 0 || particleStart[1] != 0 || particleStart[2] != 0) cout << endl << endl << endl << "###############################################" << endl << "Particle position has an issue" << endl << "###############################################" << endl;
      }

      if(verbose){
	printf("Track ipnu: %d\n",wcsimroottrack->GetIpnu());
	printf("Track parent ID: %d\n",wcsimroottrack->GetParenttype());
	printf("Track energy: %f\n", wcsimroottrack->GetE());
	printf("Track momentum: %f\n", wcsimroottrack->GetP());
	printf("Track mass: %f\n", wcsimroottrack->GetM());
      
	for (int j=0; j<3; j++){
	  printf("Track start: %d %f\n",j, wcsimroottrack->GetStart(j));
	  printf("Track dir: %d %f\n",j, wcsimroottrack->GetDir(j));
	}
      }
      }

      
    }  // End of loop over tracks
  
    double dwall, dwallBarrel, dwallCap;
    double StartRadius,StartHeight;
    if(Gamma){
      StartRadius = TMath::Sqrt(particleStop[0]*particleStop[0]+particleStop[1]*particleStop[1]);
      StartHeight = TMath::Abs(particleStop[2]);
    }
    else{
      StartRadius = TMath::Sqrt(particleStart[0]*particleStart[0]+particleStart[1]*particleStart[1]);
      StartHeight = TMath::Abs(particleStart[2]);
    }
    dwallBarrel = TankRadius - StartRadius;
    dwallCap = TankHalfHeight - StartHeight;
    dwall = dwallBarrel;//min(dwallBarrel,dwallCap);
    if(ev%1 == 0){
      cout<<"Evt #"<<ev<<endl;
      //cout<<"**************************"<<endl<<"Vertex:"<<Start_x[0]<<", "<<Start_y[0]<<", "<<Start_z[0]<<endl;
      cout<<"Dwall = "<<dwall<<", radius = "<<StartRadius<<", Height = "<<StartHeight<<endl;
      //VertexXdWall->Fill(dwall);
    }
    // Now look at the Cherenkov hits
    
    // Get the number of Cherenkov hits.
    // Note... this is *NOT* the number of photons that hit tubes.
    // It is the number of tubes hit with Cherenkov photons.
    // The number of digitized tubes will be smaller because of the threshold.
    // Each hit "raw" tube has several photon hits.  The times are recorded.
    // See chapter 5 of ../doc/DetectorDocumentation.pdf
    // for more information on the structure of the root file.
    //  
    // The following code prints out the hit times for the first 10 tubes and also
    // adds up the total pe.
    // 
    // For digitized info (one time/charge tube after a trigger) use
    // the digitized information.
    //

    int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
    int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits(); 
    int ncherenkovhits2 = 0; if(hybrid) ncherenkovhits2 = wcsimrootevent2->GetNcherenkovhits();
    int ncherenkovdigihits2 = 0;if(hybrid) ncherenkovdigihits2 = wcsimrootevent2->GetNcherenkovdigihits(); 
    
    h1->Fill(ncherenkovdigihits);
    if(verbose){
      printf("node id: %i\n", ev);
      printf("Ncherenkovhits %d\n",     ncherenkovhits);
      printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);
      printf("Ncherenkovhits2 %d\n",     ncherenkovhits2);
      printf("Ncherenkovdigihits2 %d\n", ncherenkovdigihits2);
      cout << "RAW HITS:" << endl;
    }

    if(verbose) cout << "DIGITIZED HITS:" << endl;
    //for (int index = 0 ; index < wcsimrootsuperevent->GetNumberOfEvents(); index++) 
    //{
    //wcsimrootevent = wcsimrootsuperevent->GetTrigger(index);
    //if(verbose) cout << "Sub event number = " << index << "\n";
    
    //int ncherenkovdigihits = wcsimrootevent->GetNcherenkovdigihits();
    //if(verbose) printf("Ncherenkovdigihits %d\n", ncherenkovdigihits);

    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      if(verbose) cout << "PMT Type = " << pmtType << endl;
      // Grab the big arrays of times and parent IDs
      TClonesArray *timeArray;
      if(pmtType==0) timeArray = wcsimrootevent->GetCherenkovHitTimes();
      else timeArray = wcsimrootevent2->GetCherenkovHitTimes();
      
      double particleRelativePMTpos[3];
      double totalPe = 0;
      int totalHit = 0;

      int fullHits;
      if(pmtType == 0) fullHits = ncherenkovdigihits;
      else fullHits = ncherenkovdigihits2;

      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      int unmaskedHit=0;
      for (i=0; i< fullHits ; i++)
	{
	  TObject *Hit;
	  if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	  else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);
	  WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
	    dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
	  int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId();
	  //cout<<"tube ID = "<<tubeNumber-1<<", masked = "<<maskedTube[pmtType][tubeNumber]<<endl;
	  if(maskPMT && maskedTube[pmtType][tubeNumber]){
	    //cout<<"masked"<<endl;
	    //return 0;
	    continue;
	  }
	  unmaskedHit++;
	} // End of loop over Cherenkov hits
      int nhits = unmaskedHit;
      cout<<"NHits after masking ="<<nhits<<endl;


      
      double ** essentialHitInfo_temp = new double*[nhits];
      
      // Loop through elements in the TClonesArray of WCSimRootCherenkovHits
      unmaskedHit=0;
      for (i=0; i< fullHits ; i++)
	{
	  TObject *Hit;
	  if(pmtType==0) Hit = (wcsimrootevent->GetCherenkovDigiHits())->At(i);
	  else Hit = (wcsimrootevent2->GetCherenkovDigiHits())->At(i);

	  WCSimRootCherenkovDigiHit *wcsimrootcherenkovdigihit = 
	    dynamic_cast<WCSimRootCherenkovDigiHit*>(Hit);
	  
	  int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId();
	  //cout<<"tube ID = "<<tubeNumber-1<<", masked = "<<maskedTube[pmtType][tubeNumber]<<endl;
	  if(maskPMT && maskedTube[pmtType][tubeNumber]){
	    //cout<<"masked"<<endl;
	    //return 0;
	    continue;
	  }
	  double peForTube      = wcsimrootcherenkovdigihit->GetQ();
	  double time = wcsimrootcherenkovdigihit->GetT();

	  //int tubeNumber     = wcsimrootcherenkovhit->GetTubeID();
	  //int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
	  //int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
	  WCSimRootPMT pmt;
	  if(pmtType == 0) pmt = geo->GetPMT(tubeNumber-1,false);
	  else pmt  = geo->GetPMT(tubeNumber-1,true); 
	  
	  double PMTpos[3];
	  double PMTdir[3];
	  for(int j=0;j<3;j++){
	    PMTpos[j] = pmt.GetPosition(j);
	    PMTdir[j] = pmt.GetOrientation(j);
	  }
	  
	  
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  double radius = TMath::Sqrt(PMTpos[0]*PMTpos[0] + PMTpos[1]*PMTpos[1]);
	  double phiAngle;
	  if(PMTpos[1] >= 0) phiAngle = TMath::ACos(PMTpos[0]/radius);
	  else phiAngle = TMath::Pi() + (TMath::Pi() - TMath::ACos(PMTpos[0]/radius));
	  
	  if(Gamma){//Produce the profile at the gamma conversion point, since it takes ~36 cm to convert/be visible
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStop[j];
	  }
	  else{
	    for(int j=0;j<3;j++) particleRelativePMTpos[j] = PMTpos[j] - particleStart[j];
	  }
	  
	  double vDir[3];double vOrientation[3];
	  double vDir2D[3];double vOrientation2D[3];
	  for(int j=0;j<3;j++){
	    vDir[j] = particleRelativePMTpos[j];
	    vOrientation[j] = PMTdir[j];
	    if(j<2){
	      vDir2D[j] = particleRelativePMTpos[j];
	      vOrientation2D[j] = PMTdir[j];
	    }	  
	  }
	  double Norm = TMath::Sqrt(vDir[0]*vDir[0]+vDir[1]*vDir[1]+vDir[2]*vDir[2]);
	  double tof = Norm/(lightSpeed);
	  double NormOrientation = TMath::Sqrt(vOrientation[0]*vOrientation[0]+vOrientation[1]*vOrientation[1]+vOrientation[2]*vOrientation[2]);
	  for(int j=0;j<3;j++){
	    vDir[j] /= Norm;
	    vOrientation[j] /= NormOrientation;
	  }
	  double Norm2D = TMath::Sqrt(vDir2D[0]*vDir2D[0]+vDir2D[1]*vDir2D[1]);
	  double NormOrientation2D = TMath::Sqrt(vOrientation2D[0]*vOrientation2D[0]+vOrientation2D[1]*vOrientation2D[1]);
	  for(int j=0;j<2;j++){
	    vDir2D[j] /= Norm2D;
	    vOrientation2D[j] /= NormOrientation2D;
	  }
	  double HorizontalScalarProduct = TMath::ACos(vDir2D[0]*vOrientation2D[0] + vDir2D[1]*vOrientation2D[1])*180./TMath::Pi();
	  double relativeAngle = TMath::ACos(vDir[0]*particleDir[0]+vDir[1]*particleDir[1]+vDir[2]*particleDir[2])*180./TMath::Pi();
	  //WCSimRootCherenkovHitTime HitTime = dynamic_cast<WCSimRootCherenkovHitTime>(timeArray->At(timeArrayIndex));//Assumes that the earliest time the PMT is hit is in the first element -> Everything is ordered.
	  //double time = HitTime.GetTruetime();

	  if(verbose==2){
	    if ( i < 10 ) // Only print first XX=10 tubes
	      {
		if(verbose) printf("Total pe: %1.1f times( ",peForTube);
		if(verbose) cout << ")" << endl;
		if(verbose){
		  std::cout << "Position of PMT = " << PMTpos[0] << ", " << PMTpos[1] << ", " << PMTpos[2] << endl;
		  std::cout << "Timing of the PMT hit = " << time << ", distance to vertex = " << Norm << endl;
		}
	      }
	  }
	 
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  ////////////////////////////////////////////////////////////////////////////
	  if(unmaskedHit>=fullHits){
	    cout<<"ERROR IN MASKING PMT"<<endl;
	    return 0;
	  }
	  essentialHitInfo_temp[unmaskedHit] = new double[6];
	  //essentialHitInfo_temp[unmaskedHit][0] = time+triggerTime[pmtType];//Remove the effect of trigger
	  essentialHitInfo_temp[unmaskedHit][0] = time+triggerTime[pmtType]-triggerShift[pmtType];//Remove the effect of trigger and 950 shift
	  for(int j=0;j<3;j++) essentialHitInfo_temp[unmaskedHit][j+1] = PMTpos[j];//Conversion from cm to m removed
	  essentialHitInfo_temp[unmaskedHit][4] = pmtType;
	  essentialHitInfo_temp[unmaskedHit][5] = tubeNumber;
	  //std::cout << "Relative angle = " << relativeAngle << ", time = " << time << endl;
	  
	  if(plotDigitized){
	    //plotBQ(PMTpos[3], PMTdir[3],phiAngle,relativeAngle,peForTube,pmtType);
	    //if(TMath::Abs(PMTpos[2]) < 4500){//Temp
	      if(phiAngle < TMath::Pi()/2 || phiAngle > 3*TMath::Pi()/2) ChargeProfile2D_onlyFront[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2D[pmtType]->Fill(PMTpos[1],PMTpos[2],peForTube);
	      ChargeProfile2DTop[pmtType]->Fill(PMTpos[0],PMTpos[1],peForTube);
	      ChargeProfile2DCylindricCoordinate[pmtType]->Fill(phiAngle*180/TMath::Pi(),PMTpos[2],peForTube);
	      ChargeProfile2DRelativeAngle[pmtType]->Fill(relativeAngle,PMTpos[2],peForTube);
	      
	      ChargeProfile[pmtType]->Fill(relativeAngle,peForTube);
	      HitProfile[pmtType]->Fill(relativeAngle,1);
	      
	      TimeProfile[pmtType]->Fill(time,peForTube);
	      TimeTOFProfile[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],peForTube);
	      TimeTOFProfileXTOF[pmtType]->Fill(time-tof+triggerTime[pmtType]-triggerShift[pmtType],tof,peForTube);
	      
	      ChargePerPMT[pmtType]->Fill(peForTube);


	      totalPe += peForTube;
	      totalHit ++;
	      //}
	    //if(i%50 == 0) cout << "Hit angle wrt particle direction = " << relativeAngle << ", Phi = " << phiAngle*180/TMath::Pi() << " Z position = " << PMTpos[2] << endl;
	  }
	  unmaskedHit++;
	} // End of loop over Cherenkov hits

      if(verbose) cout << "Total Pe : " << totalPe << endl;
      if(plotDigitized){
	TotalCharge[pmtType]->Fill(totalPe);
	TotalHit[pmtType]->Fill(totalHit);
	TotalChargeXdWall[pmtType]->Fill(dwall,totalPe);
      }
    
      ///////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////////////////////
      int nhitsTotal;//Number of hits in the essential hit table, that can be different than in the event because we use both PMT type hits sometimes
      if(pmtType==1){
	//Case 2: second PMT type. In that case, we will provide the list of both first and secod PMT types. Note that a copy of the first PMT type list is comtained in essentialHitIngo_pmtType0
	//Create in that case a list containing the whole hits.
	nhitsTotal=nhits_pmtType0+nhits;
	essentialHitInfo = new double*[nhitsTotal];
	//Fill it first with the type0 PMT hits.
	for (i=0; i< nhitsTotal ; i++)
	  {
	    if(i<nhits_pmtType0){
	      essentialHitInfo[i] = new double[6];
	      for(int j=0;j<6;j++) essentialHitInfo[i][j] = essentialHitInfo_pmtType0[i][j];
	      delete essentialHitInfo_pmtType0[i];
	    }
	    else{
	      essentialHitInfo[i] = new double[6];
	      for(int j=0;j<6;j++) essentialHitInfo[i][j] = essentialHitInfo_temp[i-nhits_pmtType0][j];//Remove the effect of trigger and 950 shift	      
	    }
	  }
	cout<<"Ev="<<ev<<", pmt type="<<pmtType<<",nhits pmtType0="<<nhits_pmtType0<<", nhitsTotal="<<nhitsTotal<<endl;	  
	delete essentialHitInfo_pmtType0;
      }
      else{
	nhitsTotal=nhits;
	essentialHitInfo = new double*[nhitsTotal];
	for (i=0; i< nhitsTotal ; i++)
	  {
	    essentialHitInfo[i] = new double[6];
	    for(int j=0;j<6;j++) essentialHitInfo[i][j] = essentialHitInfo_temp[i][j];//Remove the effect of trigger and 950 shift
	  }
	cout<<"Ev="<<ev<<", pmt type="<<pmtType<<",nhitsTotal="<<nhitsTotal<<endl;
	//Case 1: firt PMT type, then create the hit list as a copy of the temporary one
	if(pmtType==0){
	  essentialHitInfo_pmtType0 = new double*[nhits];
	  nhits_pmtType0=nhits;
	  for (i=0; i< nhitsTotal ; i++)
	    {
	  essentialHitInfo_pmtType0[i] = new double[6];
	      for(int j=0;j<6;j++) essentialHitInfo_pmtType0[i][j] = essentialHitInfo[i][j];//Remove the effect of trigger and 950 shift
	    }
	}
      }
      for(int i=0;i<nhits;i++){
	delete essentialHitInfo_temp[i];
      }
      delete essentialHitInfo_temp;
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




     
      for(int ibiny=1;ibiny<=timeWindowXsignal[pmtType]->GetNbinsY();ibiny++){
	for(int ibinz=1;ibinz<=timeWindowXsignal[pmtType]->GetNbinsZ();ibinz++){
	  int lowerLimitTime=timeWindowXsignal[pmtType]->GetYaxis()->GetBinLowEdge(ibiny);
	  int upperLimitTime=timeWindowXsignal[pmtType]->GetZaxis()->GetBinLowEdge(ibinz);
	  int nhitsInTime=TimeTOFProfile_currentEvent[pmtType]->Integral(TimeTOFProfile_currentEvent[pmtType]->FindBin(lowerLimitTime),TimeTOFProfile_currentEvent[pmtType]->FindBin(upperLimitTime));
	  int DR500ns = TimeTOFProfile_currentEvent[pmtType]->Integral(-600,-100);//Find dark rate in a 500ns time window out of signal.
	  double windowLength=(upperLimitTime-lowerLimitTime);
	  double averageDRInTime = DR500ns/(windowLength==0?1:windowLength);
	  double signalInTime = nhitsInTime - averageDRInTime;
	  //double optimizer = pow(signalInTime,3/2) / averagenhitsDRInTime;
	  timeWindowXsignal[pmtType]->Fill(totalHit,lowerLimitTime,upperLimitTime,signalInTime);
	  timeWindowXDR[pmtType]->Fill(totalHit,lowerLimitTime,upperLimitTime,averageDRInTime);
	}
      }

      
      //1. Read the size of the detector:
      double tankRadius = geo->GetWCCylRadius();//Conversion from cm to m removed
      double tankHeight = geo->GetWCCylLength();//Conversion from cm to m removed
      double integrationTimeWindow = 50;//in ns
      double * trueVertexPosition = new double[5];
      trueVertexPosition[0] = 0;//Since neutrino generates at t=0ns
      for(int i=0;i<3;i++){
	trueVertexPosition[i+1] = particleStart[i];//onversion to meter
      }
    

      //2. Proceed to the fit.
      double * limits = new double[4];
      if(vertexFit){
		
	
	if(nhitsTotal>0){
	  //2.a. Coarse vertex search from scratch
	  double stepSize=300;//in centimeters
	  int tolerance=60;//30
	  double ** reconstructedVertexPosition = searchVertex(tankRadius,tankHeight,integrationTimeWindow,stepSize,nhitsTotal,pmtType,trueVertexPosition,tolerance,verbose,false,stimePDFLimitsQueueNegative,stimePDFLimitsQueuePositive,false/*useDirectionality*/);
	  if(verbose==2){
	    for(int a=0;a<tolerance;a++) cout<<"Candidate vertex time = "<<reconstructedVertexPosition[a][0]<<endl;
	  }
	  ////////////////////////////////////////
#ifdef KILLFAKE
	  bool trueIncluded=false;
	  for(int a=0;a<tolerance;a++){
	    
	    double timetotrue=0;
	    double distancetotrue=0;
	    double nllDifference=0;
	    for(int j=0;j<4;j++){
	      if(j==0) timetotrue = reconstructedVertexPosition[a][j]-trueVertexPosition[j];
	      else{
		//cout<<"Distance in "<<j<<" is ="<<reconstructedVertexPositionFine3[a][j]<<","<<trueVertexPosition[j]<<endl;
		//cout<<"Distance in "<<j<<" is ="<<pow(reconstructedVertexPositionFine3[a][j]-trueVertexPosition[j],2)<<endl;
		distancetotrue+=(reconstructedVertexPosition[a][j]-trueVertexPosition[j])*(reconstructedVertexPosition[a][j]-trueVertexPosition[j]);
	      }
	    }
	    distancetotrue=TMath::Sqrt(distancetotrue);
	    if(distancetotrue < 2*stepSize && timetotrue < 2*stepSize/lightSpeed) trueIncluded=true;
	  }
	  if(!trueIncluded) continue;
#endif
	  double ** reconstructedVertexPositionFinal;
	  double stepSizeFinal=10;
	  int toleranceFinal=1;
	    if(stepByStep){
	    //2.b. Refined vertex search around candidate found previously
	    //double stepSizeFine = 300;
	    //int toleranceFine=40;
	    //for(int i=0;i<4;i++) limits[i+1] = stepSize;//particleStart[i]*1e-2;//onversion to meter
	    
	    //double ** reconstructedVertexPositionFine = searchVertexFine(reconstructedVertexPosition,limits,stepSizeFine,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFine,verbose);
	    //double ** reconstructedVertexPositionFine = searchVertexFine(reconstructedVertexPosition,limits,stepSizeFine,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFine,verbose,true,false,-100,500);
	    //if(verbose == 2){
	      //for(int a=0;a<toleranceFine;a++) cout<<"Fine candidate vertex time = "<<reconstructedVertexPositionFine[a][0]<<endl;
	      //}
	    ////////////////////////////////////////
	    
	    //2.c. More refined vertex search around candidate found previously
	    double stepSizeFine2=100;
	    int toleranceFine2=20;
	    for(int i=0;i<4;i++) limits[i] = stepSize;//particleStart[i]*1e-2;//onversion to meter
	    double ** reconstructedVertexPositionFine2 = searchVertexFine(reconstructedVertexPosition,limits,stepSizeFine2,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFine2,verbose,true,false,stimePDFLimitsQueueNegative,stimePDFLimitsQueuePositive,useDirectionality);
	    //double ** reconstructedVertexPositionFine2 = searchVertexFine(reconstructedVertexPositionFine,limits,stepSizeFine2,nhitsTotal,pmtType,trueVertexPosition,toleranceFine,toleranceFine2,1,true,false,-100,500);
	    //////////////////////////////////////////
	    
	    //2.d. Final refined vertex search around candidate found previously
	    double stepSizeFine3=20;
	    int toleranceFine3=1;
	    for(int i=0;i<4;i++) limits[i] = stepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
	    double ** reconstructedVertexPositionFine3 = searchVertexFine(reconstructedVertexPositionFine2,limits,stepSizeFine3,nhitsTotal,pmtType,trueVertexPosition,toleranceFine2,toleranceFine3,verbose,true,false,stimePDFLimitsQueueNegative,stimePDFLimitsQueuePositive,useDirectionality);
	    //double ** reconstructedVertexPositionFine3 = searchVertexFine(reconstructedVertexPositionFine2,limits,stepSizeFine3,nhitsTotal,pmtType,trueVertexPosition,toleranceFine2,toleranceFine3,verbose,true,false,-100,500);
	    //for(int i=0;i<4;i++) limits[i] = 0;//stepSizeFine;//particleStart[i]*1e-2;//onversion to meter
	    //double ** reconstructedVertexPositionFineTest = searchVertexFine(reconstructedVertexPositionFine3,limits,stepSizeFine3,nhitsTotal,pmtType,trueVertexPosition,1,toleranceFine3,2,true,false);
	    if(verbose == 1){
	      for(int a=0;a<toleranceFine3;a++) cout<<"Final candidate vertex time = "<<reconstructedVertexPositionFine3[a][0]<<", x = "<<reconstructedVertexPositionFine3[a][1]<<", y = "<<reconstructedVertexPositionFine3[a][2]<<", z="<<reconstructedVertexPositionFine3[a][3]<<endl;
	    }
	    //////////////////////////////////////////	
	    PDFnormalization_fullTimeWindow=splineIntegral(stimePDFQueue[pmtType],stimePDFLimitsQueueNegative_fullTimeWindow,stimePDFLimitsQueuePositive_fullTimeWindow);
	    for(int i=0;i<4;i++) limits[i] = 100;//stepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
	    //double ** reconstructedVertexPositionFinal = reconstructedVertexPositionFine3;//searchVertexFine(,limits,stepSizeFinal,nhitsTotal,pmtType,trueVertexPosition,toleranceFine3,toleranceFinal,verbose,true,false);
	    //double ** reconstructedVertexPositionFinal = minimizeVertex(&trueVertexPosition,limits,stepSizeFinal,nhitsTotal,pmtType,trueVertexPosition,toleranceFine3,toleranceFinal,verbose,true,false,-100,100);
	    reconstructedVertexPositionFinal = minimizeVertex(reconstructedVertexPositionFine3,limits,stepSizeFinal,nhitsTotal,pmtType,trueVertexPosition,toleranceFine3,toleranceFinal,verbose,true,false,-700,1000,useDirectionality);
	    if(verbose == 1){
	      for(int a=0;a<toleranceFinal;a++){
		cout<<"Final candidate vertex time = "<<reconstructedVertexPositionFinal[a][0]<<", x = "<<reconstructedVertexPositionFinal[a][1]<<", y = "<<reconstructedVertexPositionFinal[a][2]<<", z="<<reconstructedVertexPositionFinal[a][3]<<endl;
	      }
	    }
	    }
	  else{
	    PDFnormalization_fullTimeWindow=splineIntegral(stimePDFQueue[pmtType],stimePDFLimitsQueueNegative_fullTimeWindow,stimePDFLimitsQueuePositive_fullTimeWindow);
	    for(int i=0;i<4;i++) limits[i] = 2*stepSize;
	    reconstructedVertexPositionFinal = minimizeVertex(reconstructedVertexPosition,limits,stepSizeFinal,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFinal,verbose,true,false,-700,1000,useDirectionality);
	    //double stepSizeFine = 50;
	    //int toleranceFine=30;
	    //double ** reconstructedVertexPositionFine = minimizeVertex(reconstructedVertexPosition,limits,stepSizeFine,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFine,verbose,true,false,-100,500);
	    //for(int i=0;i<4;i++) limits[i] = stepSizeFine*4;
	    //reconstructedVertexPositionFinal = minimizeVertex(reconstructedVertexPositionFine,limits,stepSizeFinal,nhitsTotal,pmtType,trueVertexPosition,tolerance,toleranceFinal,verbose,true,false,-100,500);
	    if(verbose == 1){
	      for(int a=0;a<toleranceFinal;a++) cout<<"Final candidate vertex time = "<<reconstructedVertexPositionFinal[a][0]<<", x = "<<reconstructedVertexPositionFinal[a][1]<<", y = "<<reconstructedVertexPositionFinal[a][2]<<", z="<<reconstructedVertexPositionFinal[a][3]<<endl;
	    }
	  }

  /*
	  //Bonus:
	  //double stepSizeFinal=20;
	  //int toleranceFine3=1;
	  for(int i=0;i<4;i++) limits[i] = 0;//stepSizeFine2;//particleStart[i]*1e-2;//onversion to meter
	  double ** reconstructedVertexPositionFinal;
	  //for(int a=0;a<toleranceFine3;a++){
	  reconstructedVertexPositionFinal = searchVertexFine(reconstructedVertexPositionFine3,limits,stepSizeFine3,nhitsTotal,pmtType,trueVertexPosition,toleranceFine3,toleranceFine3,1,true,-1,2);
	  if(verbose == 1){
	    for(int a=0;a<toleranceFine3;a++){
	      //cout<<"Updated final candidate vertex time = "<<reconstructedVertexPositionFinal[a][0]<<", x = "<<reconstructedVertexPositionFinal[a][1]<<", y = "<<reconstructedVertexPositionFinal[a][2]<<", z="<<reconstructedVertexPositionFinal[a][3]<<endl;
	    }
	  }
	  //}
	//
	*/

	  //3.a. Test the NLL of the true vertex to compare with reco
	  //if(verbose) cout<<"Test of the NLL at the candidate vertex point"<<endl;
	  //for(int i=0;i<4;i++) limits[i] = 0;//stepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
	  //double ** reconstructedVertexPositionTest0 = searchVertexFine(&(reconstructedVertexPositionFine[0]),limits,stepSizeFineFine,nhitsTotal,pmtType,trueVertexPosition,1,100,verbose,true);
	
	  if(verbose) cout<<"Test of the NLL at the true vertex point"<<endl;
	  for(int i=0;i<4;i++) limits[i] = 0;//stepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
	  double ** trueVertexPositionFinal = searchVertexFine(&trueVertexPosition,limits,stepSize,nhitsTotal,pmtType,trueVertexPosition,1,1,verbose,true);
	  if(verbose) cout<<"done"<<endl;
	  //DistanceToTrueXNLL[pmtType]->Fill(distancetotrue);
	  /////////////////////////////////////////////////////////////

	  
	  // Bonus: Probe the NLL space near the true vertex:
	  if(pmtType==10 && (ev%2 == 1)){
	    for(int i=0;i<4;i++) limits[i] = 2000;//stepSizeFineFine;//particleStart[i]*1e-2;//onversion to meter
	    double stepSizeTrue=500;
	    int toleranceTrue=200;
	    cout<<"Scan"<<endl;
	    //double ** reconstructedVertexPosition = searchVertex(tankRadius,tankHeight,integrationTimeWindow,stepSize,nhitsTotal,pmtType,trueVertexPosition,tolerance,verbose,false,stimePDFLimitsQueueNegative,stimePDFLimitsQueuePositive,false/*useDirectionality*/);
	    //double ** trueVertexPositionScan = searchVertex(&trueVertexPosition,limits,stepSizeTrue,nhitsTotal,pmtType,trueVertexPosition,1,toleranceTrue,2,true,false,stimePDFLimitsQueueNegative,stimePDFLimitsQueuePositive,2/*useDirectionality*/);
	    double ** trueVertexPositionScan = searchVertexFine(&trueVertexPosition,limits,stepSizeTrue,nhitsTotal,pmtType,trueVertexPosition,1,toleranceTrue,2,true,false,-10,20,2/*useDirectionality*/);

	    for(int a=0;a<toleranceTrue;a++){
	      
	      double timetotrue=0;
	      double distancetotrue=0;
	      double nllDifference=0;
	      for(int j=0;j<4;j++){
		if(j==0) timetotrue = trueVertexPositionScan[a][j]-trueVertexPosition[j];
		else{
		  //cout<<"Distance in "<<j<<" is ="<<reconstructedVertexPositionFine3[a][j]<<","<<trueVertexPosition[j]<<endl;
		  //cout<<"Distance in "<<j<<" is ="<<pow(reconstructedVertexPositionFine3[a][j]-trueVertexPosition[j],2)<<endl;
		  distancetotrue+=(trueVertexPositionScan[a][j]-trueVertexPosition[j])*(trueVertexPositionScan[a][j]-trueVertexPosition[j]);
		}
	      }
	      distancetotrue=TMath::Sqrt(distancetotrue);
	      nllDifference=trueVertexPositionScan[a][4]-trueVertexPositionFinal[0][4];
	      fitNLLprofile_currentEvent[pmtType]->Fill(timetotrue,distancetotrue,nllDifference);
	      normfitNLLprofile_currentEvent[pmtType]->Fill(timetotrue,distancetotrue,1);
	      //distancetotrue=TMath::Sqrt(distancetotrue)*1e2;//Conversion back to cm
	    }
	  }
	  ///////////////////////////////////////////////
	
	  
	  //3.b. For final candidate, test how far it is from true vertex
	  double timetotrue=0;
	  double distancetotrue=0;
	  double nllDifference=0;

	  if(verbose) cout<<"done"<<endl;
	  for(int a=0;a<toleranceFinal;a++){
	    if(verbose) cout<<"Final candidate vertex time = "<<reconstructedVertexPositionFinal[a][0]<<endl;
	    if(a==0){
	      for(int j=0;j<4;j++){
		if(j==0) timetotrue = reconstructedVertexPositionFinal[a][j]-trueVertexPosition[j];
		else{
		  cout<<"Distance in "<<j<<" is ="<<reconstructedVertexPositionFinal[a][j]<<","<<trueVertexPosition[j]<<endl;
		  //cout<<"Distance in "<<j<<" is ="<<pow(reconstructedVertexPositionFine3[a][j]-trueVertexPosition[j],2)<<endl;
		  distancetotrue+=(reconstructedVertexPositionFinal[a][j]-trueVertexPosition[j])*(reconstructedVertexPositionFinal[a][j]-trueVertexPosition[j]);
		}
	      }
	      nllDifference=reconstructedVertexPositionFinal[a][4]-trueVertexPositionFinal[0][4];
	      //distancetotrue=TMath::Sqrt(distancetotrue)*1e2;//Conversion back to cm
	      
	    }
	  }
	  if(verbose) cout<<"done"<<endl;
	  distancetotrue=TMath::Sqrt(distancetotrue);
	  if(trueVertexPositionFinal[0][4]!=0) nllDifference/=trueVertexPositionFinal[0][4];
	  if(verbose) cout<<"Vertex found has a time difference of "<<timetotrue<<"ns and a distance of "<<distancetotrue<<"cm from the true vertex, and nll difference ="<<nllDifference<<endl;
	  if(distancetotrue>500) cout<<"code 666"<<endl;

	  fitTimeToTrue[pmtType]->Fill(timetotrue);
	  fitDistanceToTrue[pmtType]->Fill(distancetotrue);
	  fit4DdistanceXNLLdifference[pmtType]->Fill(distancetotrue,nllDifference);
	  cout << endl;
	  //////////////////////////////////////////////////////////////////
	  for(int j=0;j<4;j++){
	    /* if(j<3){
	      bsvertex[j] = reconstructedVertexPositionFinal[0][j+1];
	      mcvertex[j] = trueVertexPositionFinal[0][j+1];
	    }
	    else{
	      bsvertex[3] = reconstructedVertexPositionFinal[0][0];
	      mcvertex[3] = trueVertexPositionFinal[0][0];
	    }*/
	    if(j==0){
	      bsvertex[0] = reconstructedVertexPositionFinal[0][3];
	      mcvertex[0] = trueVertexPosition[3];
	    }
	    else if(j==3){
	      bsvertex[3] = reconstructedVertexPositionFinal[0][0];
	      mcvertex[3] = trueVertexPosition[0];
	    }
	    else{
	      bsvertex[j] = reconstructedVertexPositionFinal[0][j];
	      mcvertex[j] = trueVertexPosition[j];
	    }
	    
	  }
	  for(int j=0;j<3;j++){
	    mcdir[j]=particleDir[j];
	  }
	  bsnhit[0] = nhitsTotal;
	  pmtTypeFill = pmtType;

	  bsnhit_intime[pmtType]=0;
	  for(int ihit = 0; ihit < nhitsTotal; ihit++){
	    double * info = essentialHitInfo[ihit];
	    double hitTime = info[0];
	    double hitPosition[3];
	    double distance = 0;
	    for(int i=0;i<3;i++){
	      hitPosition[i] = info[i+1] - reconstructedVertexPositionFinal[0][i+1];
	      distance+=pow(hitPosition[i],2);
	    }
	    distance=TMath::Sqrt(distance);
	    double tof = distance / lightSpeed;
	    double residual = hitTime - tof - reconstructedVertexPositionFinal[0][0];
	    if(residual > stimePDFLimitsQueueNegative && residual < stimePDFLimitsQueuePositive){
	      bsnhit_intime[pmtType]++;
	    }
	  }
	  bstree->Fill();
	}
      }
      cout<<"debug, wil delete"<<endl;
      for(int i=0;i<nhitsTotal;i++){
	delete essentialHitInfo[i];
      }
      delete trueVertexPosition;
      delete essentialHitInfo;
      delete limits;
      cout<<"debug, go to next"<<endl;
    }
    // Get the number of digitized hits
    // Loop over sub events
    
    // reinitialize super event between loops.
    wcsimrootsuperevent->ReInitialize();
    if(hybrid) wcsimrootsuperevent2->ReInitialize();

#ifdef DRAW
    for(int i=0;i<nPMTtypes;i++){
      TimeTOFProfile_currentEvent[i]->Reset();
      if(ev-startEvent<10){
	fitNLLprofile_currentEvent[i]->Divide(normfitNLLprofile_currentEvent[i]);
	fitNLLprofile_currentEvent[i]->Write(Form("fitNLLprofile_event%d_pmtType%d",ev,i));
      }
      fitNLLprofile_currentEvent[i]->Reset();
      normfitNLLprofile_currentEvent[i]->Reset();
    }
#endif
  } // End of loop over events
  //  TCanvas c1("c1"); 
  /*
  float win_scale = 0.75;
  int n_wide(2);
  int n_high(2);
  TCanvas* c1 = new TCanvas("c1", "First canvas", 500*n_wide*win_scale, 500*n_high*win_scale);
  c1->Draw();
  c1->Divide(2,2);
  c1->cd(1); hvtx0->Draw();
  c1->cd(2); hvtx1->Draw();
  c1->cd(3); hvtx2->Draw();
  c1->cd(4); h1->Draw();
  */
  std::cout<<"num_trig "<<num_trig<<"\n";
  if(vertexFit) bstree->Write();
  else{
    for(int i=0;i<nPMTtypes;i++){

      ChargeProfile2D_onlyFront[i]->Write();
      ChargeProfile2D[i]->Write();
      ChargeProfile2DTop[i]->Write();
      ChargeProfile2DCylindricCoordinate[i]->Write();
      ChargeProfile2DRelativeAngle[i]->Write();
      ChargeProfile[i]->Write();
      HitProfile[i]->Write();
      TimeProfile[i]->Write();
      TimeTOFProfile[i]->Write();
      TimeTOFProfileXTOF[i]->Write();
      ChargePerPMT[i]->Write();
      TotalCharge[i]->Write();
      TotalHit[i]->Write();
      TotalChargeXdWall[i]->Write();//Fill(dWall,totalPe);
      fitTimeToTrue[i]->Write();
      fitDistanceToTrue[i]->Write();
      fit4DdistanceXNLLdifference[i]->Write();
      for(int ibinx=1;ibinx<=timeWindowXsignal[i]->GetNbinsX();ibinx++){
	for(int ibiny=1;ibiny<=timeWindowXsignal[i]->GetNbinsY();ibiny++){
	  for(int ibinz=1;ibinz<=timeWindowXsignal[i]->GetNbinsZ();ibinz++){
	    double signal = timeWindowXsignal[i]->GetBinContent(ibinx,ibiny,ibinz);
	    double DR = timeWindowXDR[i]->GetBinContent(ibinx,ibiny,ibinz);
	    double optimizer = pow(signal,3/2);
	    if(DR!=0) optimizer/=DR;
	    timeWindowXsignal[i]->SetBinContent(ibinx,ibiny,ibinz,optimizer);
	  }
	}
      }
      timeWindowXsignal[i]->Write(); 
    }
  }
  
  outfile->Close();
  
  return 0;
 }
