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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMinuit.h>
#include <TFitter.h>
using namespace std;
const int dimNoise=3;
const int dimSignal=3;
const int ntuple=5;
vector < array < double,ntuple > > dataID;//use array of c++ 11 as we cannot create a vector of traditional arrays e.g. double[5].
double CredibleInterval=0.68;

double SignalToNoise(double PC, double DR, double E, double n[dimNoise], double s[dimSignal]){//Energy should be in MeV and DR in kHz
  double PEperMeV=8;//For the 20% PC case.
  double signalHits=PEperMeV*E*(PC/20.);
  double noiseHits=250*(PC/20.)*(DR/8.4);//PC*1e3*DR*1e3/(1.5e-6);//For the 20% PC case, 8.4 kHz.
  //double signal_to_noise=noiseHits/(signal==0?1:pow(signal,3./2.));
  double noise, signal = 0;
  for(int i=0;i<dimNoise;i++){
    noise+=n[i]*pow(TMath::Sqrt(noiseHits),i);
  }
  for(int i=0;i<dimSignal;i++){
    signal+=s[i]*pow(TMath::Sqrt(signalHits),i);
  //signal+=s[i]*pow(signalHits,i);
  }
  double signal_to_noise=0;
  signal_to_noise=noise;
  if(signal!=0) signal_to_noise/=signal;

  //double signal_to_noise_sum=0;
  //for(int i=0;i<dimSignal;i++){
  //signal_to_noise_sum+=s[i]*pow(signal_to_noise,i);
  //}

  //pow(signal,3./2.)/noiseHits;//signal/(noiseHits==0?1:TMath::Sqrt(noiseHits));
  //cout<<"Energy ="<<E<<", signal/sqrt(noise)="<<signal_to_noise<<endl;
  return signal_to_noise;
}
/*
//Assumes that function is of the form (a + b*sqrt(N) + c*sqrt(N)^2)/S
double calculateRMS(double a0, double a1, double a2, int pmtType){
  double rms_vertex=0;
  //int nbins_signal_to_noise=10;
  //double max_signal_to_noise=-1e9;
  //double min_signal_to_noise=1e9;
  double binWidth=20;//cm;
  int nbins_vRes=0;
  double max_vRes=-1e9;
  double min_vRes=1e9;

  cout<<"Let's calculate the RMS for parameters="<<a0<<", "<<a1<<", "<<a2<<endl;  
  //First, determine the higher signal_to_noise value.
  for(int i=0;i<dataID.size();i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double pmt=data[3];
    if(pmt!=pmtType) continue;
    double vRes=data[4];
    double signal_to_noise=SignalToNoise(pc,dr,e,a0,a1,a2);
    if(i%10==0) cout<<"STN="<<signal_to_noise<<endl;
    //if(signal_to_noise>max_signal_to_noise) max_signal_to_noise=signal_to_noise;
    //if(signal_to_noise<min_signal_to_noise) min_signal_to_noise=signal_to_noise;
    if(vRes>max_vRes) max_vRes=vRes;
    if(vRes<min_vRes) min_vRes=vRes;
  }

  //Based on this highest value, adapt the binning
  //double binWidth=(max_signal_to_noise*1.01-min_signal_to_noise)/nbins_signal_to_noise;
  //cout<<"bin width = "<<binWidth<<", min = "<<min_signal_to_noise<<", max = "<<max_signal_to_noise<<endl;
  nbins_vRes=(max_vRes*1.01-min_vRes)/binWidth;
  //double binWidth=(max_vRes*1.01-min_vRes)/nbins_vRes;
  cout<<"bin width = "<<binWidth<<", min = "<<min_vRes<<", max = "<<max_vRes<<", number of bins ="<<nbins_vRes<<endl;

  //Then, create a collection of the vertex resolution point in each bin.
  //Also, calculate the mean of each bin.
  //double mean_vRes_binned[nbins_signal_to_noise]={0.};
  //int nPoints[nbins_signal_to_noise]={0};
  double mean_vRes_binned[nbins_vRes]={0.};
  int nPoints[nbins_vRes]={0};
  for(int i=0;i<dataID.size();i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double pmt=data[3];
    if(pmt!=pmtType) continue;
    double vRes=data[4];
    double signal_to_noise=SignalToNoise(pc,dr,e,a0,a1,a2);
    //int binID=((int) ((signal_to_noise-min_signal_to_noise) / binWidth));
    //cout<<"bin ID ="<<binID<<", STN = "<<signal_to_noise<<", vres = "<<vRes<<endl;
    //mean_vRes_binned[binID]+=vRes;
    //nPoints[binID]++;
    //mean_vRes_binned[binID]+=vRes;
    //nPoints[binID]++;
    int binID=((int) ((vRes-min_vRes) / binWidth));
    cout<<"bin ID ="<<binID<<", STN = "<<signal_to_noise<<", vres = "<<vRes<<endl;
    mean_vRes_binned[binID]+=signal_to_noise;
    nPoints[binID]++;
  }
  //Estimate the mean:
  //for(int i=0;i<nbins_signal_to_noise;i++){
  for(int i=0;i<nbins_vRes;i++){
    if(nPoints[i]!=0) mean_vRes_binned[i]/=nPoints[i];
    cout<<"mean sqrt(N)/S= "<<mean_vRes_binned[i]<<endl;
  }
  
  //Then calculate the RMS:
  for(int i=0;i<dataID.size();i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double pmt=data[3];
    if(pmt!=pmtType) continue;
    double vRes=data[4];
    double signal_to_noise=SignalToNoise(pc,dr,e,a0,a1,a2);
    //int binID=((int) ((signal_to_noise - min_signal_to_noise) / binWidth));
    //cout<<"bin ID ="<<binID<<", STN = "<<signal_to_noise<<", vertex res = "<<vRes<<", diff res = "<<vRes - mean_vRes_binned[binID]<<endl;
    //rms_vertex+=pow( (vRes - mean_vRes_binned[binID] ),2);
    int binID=((int) ((vRes - min_vRes) / binWidth));
    cout<<"bin ID ="<<binID<<", STN = "<<signal_to_noise<<", vertex res = "<<vRes<<", diff res = "<<vRes - mean_vRes_binned[binID]<<endl;
    rms_vertex+=pow( (signal_to_noise - mean_vRes_binned[binID] ),2);
    //rms_vRes_binned[binID]+=pow( (vRes-mean_vRes_binned[binID] ),2);
  }
  int npoints_total=dataID.size();
  rms_vertex=TMath::Sqrt(rms_vertex);
  cout<<"RMS before division = "<<rms_vertex;  
  if(npoints_total!=0) rms_vertex/=npoints_total;
  cout<<", RMS = "<<rms_vertex<<", npoints = "<<npoints_total<<endl;
  
  return rms_vertex;
}
*/

//Assumes that function is of the form (a + b*sqrt(N) + c*sqrt(N)^2)/S
double calculateRMS(double n[dimNoise], double s[dimSignal], int pmtType, bool drawOption=false){

  //The idea is to:
  //1. Sort the points by sqrt(N)/S.
  //2. Make average over groups of successive 4 hits.
  //3. Spline these points.
  //4. Calculate the RMS by checking, for each point, the distance wrt this points.

  //Other idea: assume that the law should be linear between a+bsqrt(N) + c sqrt(N) / S.
  //In that case, fit the points with a straight line.
  //To do so, the issue is: which error to take for the TGraphErrors. Take an error of +-1cm for vertex
  TGraphErrors * g;
  int npoints=0;
  int size = (int) dataID.size();
  for(int i=0;i<size;i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double pmt=data[3];
    double vertex_resolution=data[4];
    if(pmt!=pmtType) continue;
    //if(vertex_resolution<100 || vertex_resolution>300) continue;
    if(vertex_resolution>300) continue;
    //if(vertex_resolution<100) continue;
    //if(e>5) continue;
    npoints++;
  }
  double * x = new double[npoints];//={0.};
  double * ex = new double[npoints];//={0.};
  double * y = new double[npoints];//={0.};
  double * ey = new double[npoints];//={0.};
  double max_signal_to_noise=-1e9;
  double min_signal_to_noise=1e9;

  int count=0;
  for(int i=0;i<dataID.size();i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double pmt=data[3];
    double vertex_resolution=data[4];
    if(pmt!=pmtType) continue;
    //if(vertex_resolution<100 || vertex_resolution>300) continue;
    if(vertex_resolution>300) continue;
    //if(e>5) continue;
    //if(vertex_resolution>
    //double signal_to_noise=SignalToNoise(pc,dr,e,a0,a1,a2,a3);
    double signal_to_noise=SignalToNoise(pc,dr,e,n,s);
    x[count]=signal_to_noise;
    y[count]=vertex_resolution;
    ex[count]=0;
    ey[count]=10;
    if(signal_to_noise>max_signal_to_noise) max_signal_to_noise=signal_to_noise;
    if(signal_to_noise<min_signal_to_noise) min_signal_to_noise=signal_to_noise;
    count++;
  }
  g = new TGraphErrors(npoints,x,y,ex,ey);

  //Then, we should fit it with a straight line.
  TF1 * f = new TF1("f","pol2",min_signal_to_noise,max_signal_to_noise);
  if(!drawOption) g->Fit("f","Q0");
  else{
    TCanvas * c = new TCanvas();
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+2);
    g->GetXaxis()->SetTitle("");
    g->GetYaxis()->SetTitle("Vertex resolution (cm)");
    //g->GetXaxis()->SetRangeUser(3,15);
    g->GetYaxis()->SetRangeUser(0,300);
    g->Draw("AP");
    g->Fit("f","");
    f->Draw("same");
    c->SaveAs("plots/lawFit.eps");
  }

  //2. Now, just, provide the chi2 of the fit.
  double chi2=f->GetChisquare();
  double chi2_red=chi2/npoints;

  double rms_vertex=0;
  for(int i=0;i<npoints;i++){
    double diff=y[i]-f->Eval(x[i]);
    rms_vertex+=diff*diff;
  }
  rms_vertex=TMath::Sqrt(rms_vertex)/npoints;
  //cout<<"For points: "<<a0<<", "<<a1<<", "<<a2<<", "<<a3<<endl;  
  cout<<"For points: noise: ";
  for(int i=0;i<dimNoise;i++) cout<<n[i]<<", ";
  cout<<" signal: ";
  for(int i=0;i<dimSignal;i++) cout<<s[i]<<", ";
  cout<<endl;
  cout<<"Chi2 = "<<chi2<<", reduced = "<<chi2_red<<", rms = "<<rms_vertex<<endl<<endl;
  delete x, y, ex, ey;

  return chi2_red;//rms_vertex;
  //return rms_vertex;
}

//Look that the RMS of the function is as small as possible
void minuitLikelihood(int& nDim, double * gout, double & RMS, double par[], int flg){
  double noise[dimNoise]={0.};
  double signal[dimSignal]={0.};
  for(int i=0;i<dimNoise;i++) noise[i]=par[i];
  for(int i=0;i<dimSignal;i++) signal[i]=par[i+dimNoise];
  RMS=calculateRMS(noise,signal,par[dimNoise+dimSignal]);//
}



//here is the function where minimization gonna take place. The loop is inside this function
void minimizeNoise(int pmtType, double n[dimNoise], double s[dimSignal]){
  
  TFitter * minimizer = new TFitter(dimNoise+dimSignal+1);//4=nb de params?
  TMinuit * minuit = minimizer->GetMinuit();
  double arglist[20];
  int err=0;
  double p1=1;
  minimizer->ExecuteCommand("SET PRINTOUT",&p1,1);//quiet mode
  minuit->SetErrorDef(1);
  arglist[0]=2;
  minuit->mnexcm("SET STR",arglist,1,err);//set strategy sets the number of derivative estimated. in the present case (2), this is the highest number of it with a counter-balance: very slow!
  minimizer->SetFCN(minuitLikelihood);//here is function to minimize

  double Mig[2]={3e4,1e-6};//maxcalls and tolerance
  double Imp2[1]={10000};
  for(int i=0;i<dimNoise;i++) minimizer->SetParameter(i,Form("n%d",i),i==1?1.:0.,1e-3,0.,1e3);
  for(int i=0;i<dimSignal;i++) minimizer->SetParameter(i+dimNoise,Form("s%d",i+dimNoise),i==2?1.:0.,1e-3,-1e3,1e3);
  minimizer->SetParameter(dimNoise+dimSignal+1,"pmtType",pmtType,pmtType,pmtType-1,pmtType+1);
  minimizer->FixParameter(dimNoise+dimSignal+1);

  for(int i=0;i<dimSignal;i++) minimizer->FixParameter(dimNoise+i);
  Mig[1]=1e-8;minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e-1;minimizer->ExecuteCommand("MIGRAD",Mig,2);
  /*
  for(int i=0;i<dimSignal;i++) minimizer->ReleaseParameter(i+dimNoise);
  Mig[1]=1e-2;minimizer->ExecuteCommand("MIGRAD",Mig,2);

  
  for(int i=0;i<dimNoise;i++) minimizer->FixParameter(i);
  for(int i=0;i<dimSignal;i++) minimizer->ReleaseParameter(i+dimNoise);
  Mig[1]=1e-6;minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e-2;minimizer->ExecuteCommand("MIGRAD",Mig,2);
  
  for(int i=0;i<dimNoise;i++) minimizer->ReleaseParameter(i);
  Mig[1]=1e-6;minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e-2;minimizer->ExecuteCommand("MIGRAD",Mig,2);
  
  for(int i=0;i<dimNoise;i++) minimizer->FixParameter(i);
  minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e0;//tolerance
  minimizer->ExecuteCommand("MIGRAD",Mig,2);
  
  for(int i=0;i<dimSignal;i++) minimizer->FixParameter(dimNoise+i);
  for(int i=0;i<dimNoise;i++) minimizer->ReleaseParameter(i);
  minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e0;minimizer->ExecuteCommand("MIGRAD",Mig,2);

  for(int i=0;i<dimSignal;i++) minimizer->ReleaseParameter(i+dimNoise);
  minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e0;minimizer->ExecuteCommand("MIGRAD",Mig,2);
  */
  // for(int i=0;i<dimNoise;i++) minimizer->FixParameter(i);
  //for(int i=0;i<dimSignal;i++) minimizer->ReleaseParameter(dimNoise+i);
  //minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  //Mig[1]=1e0;minimizer->ExecuteCommand("MIGRAD",Mig,2);

  //for(int i=0;i<dimNoise;i++) minimizer->ReleaseParameter(i);
  //minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  //Mig[1]=1e0;minimizer->ExecuteCommand("MIGRAD",Mig,2);
  
  //double a * = double[3];
  for(int i=0;i<dimNoise;i++) n[i] = minimizer->GetParameter(i);
  for(int i=0;i<dimSignal;i++) s[i] = minimizer->GetParameter(i+dimNoise);
  cout<<"Minimization finished"<<endl;
}


bool migrateFromOVtoFV(double dwall, double dwall_reconstructed, double dwall_FV_definition){
  //cout<<"FV cut is set at a distance = "<<dwall_FV_definition<<", Event dwall = "<<dwall<<", reconstructed dwall = "<<dwall_reconstructed<<endl;
  if(dwall<dwall_FV_definition){//event which are truly out of FV
    if(dwall_reconstructed>=dwall_FV_definition) return true;
  }
  return false;
}

double vRes(TH1D * hvertex, double CI = CredibleInterval){
  double integralSoFar = 0;
  double integralTotal = hvertex->Integral();
  double integral = 0;
  int ibin1D=1;
  while(integral<CI && ibin1D<1e4){
    integralSoFar+=hvertex->GetBinContent(ibin1D);
    integral=integralSoFar/(integralTotal==0?1:integralTotal);
    ibin1D++;
  }
  return hvertex->GetBinCenter(ibin1D);
}

double vResAndError(TH1D * hvertex, double CI = CredibleInterval){
  double integralSoFar = 0;
  double integralTotal = hvertex->Integral();
  if(integralTotal==0) return 0;
  double integral = 0;
  int ibin1D=1;
  double error=0;
  while(integral<CI && ibin1D<1e4){
    integralSoFar+=hvertex->GetBinContent(ibin1D);
    integral=integralSoFar/(integralTotal==0?1:integralTotal);
    ibin1D++;
  }
  double res=hvertex->GetBinCenter(ibin1D);
  //cout<<"res="<<res<<", "<<integralSoFar<<", "<<integralTotal<<endl;
  error=integral*TMath::Sqrt(pow(1/integralSoFar,2)+pow(1/integralTotal,2));
  //cout<<"Error on integral = "<<error<<endl;
  double error_vertex_inf=vRes(hvertex,CI-error);
  double error_vertex_sup=vRes(hvertex,CI+error);
  double error_full=TMath::Max(TMath::Abs(res-error_vertex_inf),TMath::Abs(res-error_vertex_sup));
  //cout<<"Error on position inf = "<<error_vertex_inf<<", "<<error_vertex_sup<<", full = "<<error_full<<endl;
  error_full+=hvertex->GetBinWidth(ibin1D);//Add bin width
  return error_full;
}

int vertexResolution(double energy, double * vertex_resolution, double * efficiency, double * misIDFV, const int nPMTtypes,bool monoEnergy=true, double photocov=40., double darkrate=4.2){

   //Open the file

  //BIG PMT
  //TFile *file = new TFile("hybrid_e10.0MeV_0.0_PC20+5per_center.root");
  //TFile *file = new TFile("hybrid_e10.0MeV_0.0_PC25per.root");
  //TFile *file = new TFile("hybrid_e10.0MeV_8.4_PC25per.root");
  //TFile *file = new TFile("hybrid_e10.0MeV_8.4_PC25per_oldPDF.root");
  //TFile *file = new TFile("hybrid_e10.0MeV_8.4_PC25per_newPDF.root");
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4_PC25per.root",energy));
  //TFile *file = new TFile(Form("hybrid_e10.0MeV_8.4_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e10.0MeV_8.4_PC25per.root",energy));
  //TFile *file = new TFile(Form("hybrid_e10.0MeV_8.4_PC25per_balonly.root",energy));
  //TFile *file = new TFile(Form("hybrid_e10.0MeV_8.4_PC25per_balonly_triggeradded.root",energy));
  //TFile *file = new TFile(Form("hybrid_e4.0MeV_8.4_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e4.0MeV_8.4_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("test.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4_PC40per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4kHzAnd50Hz_PC25per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_8.4_PC30per_together.root",energy));
  //TFile *file = new TFile(Form("hybrid_e%2.1fMeV_9.3_PC25per_together.root",energy));
  //TFile *file = new TFile("hybrid_e10.0MeV_0.0_PC25per.root");
  //TFile *file = new TFile("hybrid_e3.0MeV_8.4_PC25per.root");

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged3.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_minimizerTimeAlso3.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_minimizerTimeAlso4_startAtRec.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_minimizerTimeAlso4_startAtRec_fullPDF.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTBLPDF.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_DRscaling.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_noDRscaling_highStat.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_noDRscaling_highStat.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_minimizerTimeAlso4_startAtRec_fullPDF_noScaling.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_minimizerTimeAlso4_startAtTrue.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_new_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_directionality6_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_removeFails.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_removeFails_directionalityAll.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityCut.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_removeFails_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_removeFails_directionalityAll_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF3MeV_faster_removeFails_directionalityAll_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDFcontinue3_faster_removeFails_directionalityAll_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF10MeV_faster_removeFails_directionalityAll_10kmPMT.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDFcontinue_faster_removeFails_directionalityAll.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_signalDir0.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_new.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_new.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTBLPDF_new.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_new_50Hz.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTsignalOnlyPDF_new_50Hz.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMT3MeVPDF_new_50Hz.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_largerHitWindow.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_directionality_new2.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTBLPDF_allMIGRAD.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTBLPDF_2msearch.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_averaged50cm.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-1.5and3nsTimeWindow_merged.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_new2.root",((int) energy) ));






  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_40blonly.root",((int) energy) ));

  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_halfDR.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_40blonly_halfDR.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_40blonly_halfMasked.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_10kmPMT_halfMasked.root",((int) energy) ));
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_halfDRonlyforBL.root",((int) energy) ));//20% B&L 4.2 kHz + 5k mPMT 100z
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_10kmPMT.root",((int) energy) ));//20% B&L 8.4 kHz + 10k mPMT 100z


  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_0.8masked.root",((int) energy) ));//20% B&L 4.2 kHz + 10k mPMT 100z
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_0.4masked.root",((int) energy) ));//20% B&L 4.2 kHz + 10k mPMT 100z
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_10kmPMT_halfDR.root",((int) energy) ));//20% B&L 4.2 kHz + 10k mPMT 50z

  
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_%dMeV_-3and4nsTimeWindow_merged_fullPDF_mPMTPDF_faster_directionalityAll_new2_10kmPMT_test.root",((int) energy) ));
  //TFile *file = new TFile("hybrid_e10.0MeV_0.0_PC25per_parameterCardv1.root");
  //TFile *file = new TFile(Form("wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_fulltank_0hitstrigger_10000_rec_nomask.root",((int) energy) ));//20% B&L 8.4 kHz + 10k mPMT 100z
  //TFile *file = new TFile(Form("output_wcsim_Ben.root",((int) energy) ));//20% B&L 8.4 kHz + 10k mPMT 100z
  //TFile *file = new TFile(Form("benNano_HKHybrid_20PC_10MeV_-3and4nsTimeWindow_tmp.root",((int) energy) ));//20% B&L 8.4 kHz + 10k mPMT 100z
  //TFile *file = new TFile(Form("output_wcsim_IzumiyamaSan.root",((int) energy) ));//20% B&L 8.4 kHz + 10k mPMT 100z
  //cout<<photocov<<", DR="<<darkrate<<endl;
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_fulltank_directionality_10MeV_test.root",((int) energy) ));


  //Official: for studies close to the wall.
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_directionality7_ID10_%dMeV.root",((int) energy) ));//Case with DR using plateau of the distributions, TGraph2D, and new 1/r correction from TF1
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_fulltank_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_fulltank_0.75PC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_closewall_0.75PC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk20pc_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_closewall_0.5PC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_fulltank_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hk40pc_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_0mPMTPC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_0.3mPMTPC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_0.5mPMTPC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_fulltank_0.3mPMTPC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_fulltank_0.5mPMTPC_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_fulltank_nodirectionality_%dMeV.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAF_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));
 //TFile *file = new TFile(Form("LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_nodirectionality_%dMeV.root",((int) energy) ));

  //LEAFv2
  //TFile *file = new TFile(Form("LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzbl_%dMeV_tmp0.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_nodirectionality_%dMeV_v2.root",((int) energy) ));

  //Official plot for big PMTs:
  TFile *file = new TFile(Form("LEAF_hk%.0fpc_%.1fkHzbl_%dMeV.root",photocov,darkrate,((int) energy) ));

  //cout<<"We will try to open file = "<<file->GetName()<<endl;
  //TFile *file = new TFile(Form("LEAF_hk40pc_6.0kHzbl_%dMeV.root",((int) energy) ));

  file->cd();
  double tankRadius = 7080/2;//In cm
  double tankHeight = 5480/2;//In cm
  bool verbose=true;
  bool bug=true;//false;//true;
  double radioactive_dwall=100;//in cm, how much the radioactive background spread in the detector
  //Get the a pointer to the tree from the file
  TTree *bstree = (TTree*)file->Get("Reduced");

  //Variables declaration
  double lf_time; double lf_vertex[3];
  double lf_dir[3];
  /*
  std::vector<float> true_vertex_X;
  std::vector<float> true_vertex_Y;
  std::vector<float> true_vertex_Z;
  std::vector<float> true_vertex_T;
  std::vector<float> true_dir_X;
  std::vector<float> true_dir_Y;
  std::vector<float> true_dir_Z;
  */
  std::vector<float> *true_vertex_X = 0;
  std::vector<float> *true_vertex_Y = 0;
  std::vector<float> *true_vertex_Z = 0;
  std::vector<float> *true_vertex_T = 0;
  std::vector<float> *true_dir_X = 0;
  std::vector<float> *true_dir_Y = 0;
  std::vector<float> *true_dir_Z = 0;
  
  double lf_good;
  double lf_dirKS;
  int lf_intime;
  int ID_hits_200;
  int ID_hits;
  int pmtT=0;
  int eventId;
  int triggerId;
  
  //15 MeV 10 MeV 8 MeV
  int nbin_vertex = 6000;//Number of bins used for vertex resolution
  double max_vertex = 6000;
  int nbin_direction = 100;//Number of bins used for angle resolution
  double max_direction = 180;
  int nbin_nhits = 400;
  double max_nhits = 400;
  int nbin_time = 400;
  double max_time = 50;
  int nbin_distance = 10;
  double min_distance = 0.;
  double max_distance = 200;
  double nbin_tanksize = 760;
  double min_tanksize = -3800;
  double max_tanksize = 3800;
  int nbin_goodness = 100;//Number of bins used for vertex resolution
  double max_goodness = 2;
  
  //Histo declaration
  TH2D * hvertexXnhits_intime[nPMTtypes];
  TH2D * hvertexXgoodness[nPMTtypes];
  TH2D * hvertexXnhits[nPMTtypes];
  TH2D * hvertexXdwall[nPMTtypes];
  TH2D * hvertexXtowall[nPMTtypes];
  TH2D * hdwallXdwallrec[nPMTtypes];
 
  TH2D * hdwallXnhits[nPMTtypes];
  TH2D * htowallXnhits[nPMTtypes];
  
  TH2D * hvertexLongitudinalXTime[nPMTtypes];
  TH1D * hvertexLongitudinal[nPMTtypes];
  TH1D * hvertexOrtho[nPMTtypes];
  TH1D * hvertex[nPMTtypes];
  TH1D * hradius[nPMTtypes];
  TH1D * hdirection[nPMTtypes];
  TH1D * hnhits[nPMTtypes];
  double nVertexFound[nPMTtypes];
  double nEvents[nPMTtypes];
  double nEventsTotal[nPMTtypes];
  TH2D * hvertex2D[nPMTtypes];
  TH2D * htruevertex2D[nPMTtypes];
  
  TH1D * hvertexMisIDFV[nPMTtypes]; TH1D * hvertexMisIDFV_total[nPMTtypes];
  double nMisIDFV[nPMTtypes]; double nMisIDFV_total[nPMTtypes];

  TH1D * hvertexMisIDExtCut[nPMTtypes]; TH1D * hvertexMisIDExtCut_total[nPMTtypes];
  double nMisIDExtCut[nPMTtypes]; double nMisIDExtCut_total[nPMTtypes];

  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    nVertexFound[pmtType]=0;
    nEventsTotal[pmtType]=0;
    nEvents[pmtType]=0;

    nMisIDFV[pmtType]=0;
    nMisIDFV_total[pmtType]=0;

    nMisIDExtCut[pmtType]=0;
    nMisIDExtCut_total[pmtType]=0;


    hvertexXnhits_intime[pmtType] = new TH2D(Form("hvertexXnhits_intime_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_nhits,0.,max_nhits);
    hvertexXnhits_intime[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexXnhits_intime[pmtType]->GetYaxis()->SetTitle("number of hits");

    hvertexXgoodness[pmtType] = new TH2D(Form("hvertexXgoddness_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_goodness,0.,max_goodness);
    hvertexXgoodness[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexXgoodness[pmtType]->GetYaxis()->SetTitle("goodness of fit");

    hvertexXnhits[pmtType] = new TH2D(Form("hvertexXnhits_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_nhits,0,max_nhits);
    hvertexXnhits[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexXnhits[pmtType]->GetYaxis()->SetTitle("number of hits");

    hdwallXnhits[pmtType] = new TH2D(Form("hdwallXnhits_pmtType%d",pmtType),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
    hdwallXnhits[pmtType]->GetXaxis()->SetTitle("dwall (cm)");
    hdwallXnhits[pmtType]->GetYaxis()->SetTitle("number of hits");

    htowallXnhits[pmtType] = new TH2D(Form("htowallXnhits_pmtType%d",pmtType),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
    htowallXnhits[pmtType]->GetXaxis()->SetTitle("towall (cm)");
    htowallXnhits[pmtType]->GetYaxis()->SetTitle("number of hits");

   
    hdwallXdwallrec[pmtType] = new TH2D(Form("hdwallXdwallrec_pmtType%d",pmtType),"",nbin_distance,0,max_distance,20,-1,1);
    hdwallXdwallrec[pmtType]->GetXaxis()->SetTitle("dwall true (cm)");
    hdwallXdwallrec[pmtType]->GetYaxis()->SetTitle("dwall rec  - dwall true / dwall true");

    //hdwallXpmt_to_vertex[pmtType] = new TH2D(Form("hdwallXpmt_to_vertex_pmtType%d",pmtType),"",nbin_distance,0,max_distance,20,-1,1);
    //hdwallXpmt_to_vertex[pmtType]->GetXaxis()->SetTitle("dwall true (cm)");
    //hdwallXpmt_to_vertex[pmtType]->GetYaxis()->SetTitle("dwall rec  - dwall true / dwall true");

    hvertexXdwall[pmtType] = new TH2D(Form("hvertexXdwall_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
    hvertexXdwall[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexXdwall[pmtType]->GetYaxis()->SetTitle("dwall (cm)");

    hvertexXtowall[pmtType] = new TH2D(Form("hvertexXtowall_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
    hvertexXtowall[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexXtowall[pmtType]->GetYaxis()->SetTitle("towall (cm)");

    hvertexLongitudinalXTime[pmtType] = new TH2D(Form("hvertexLongitudinalXTime_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex,nbin_time,-max_time,max_time);
    hvertexLongitudinalXTime[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertexLongitudinalXTime[pmtType]->GetYaxis()->SetTitle("dwall (cm)");
			     
    hvertexLongitudinal[pmtType] = new TH1D(Form("hvertexLongitudinal_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex);
    hvertexLongitudinal[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    
    hvertexOrtho[pmtType] = new TH1D(Form("hvertexOrtho_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex);
    hvertexOrtho[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hvertex[pmtType] = new TH1D(Form("hvertex_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex);
    hvertex[pmtType]->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
    hradius[pmtType] = new TH1D(Form("hradius_pmtType%d",pmtType),"",nbin_vertex,0,max_vertex*max_vertex);
    hradius[pmtType]->GetXaxis()->SetTitle("R2 (cm2)");
    hdirection[pmtType] = new TH1D(Form("hdirection_pmtType%d",pmtType),"",nbin_direction,0,max_direction);
    hdirection[pmtType]->GetXaxis()->SetTitle("angle from rec. to mc. direction (#circ)");
    hnhits[pmtType] = new TH1D(Form("hnhits_pmtType%d",pmtType),"",nbin_nhits,0,max_nhits);
    hnhits[pmtType]->GetXaxis()->SetTitle("number of hits");

    hvertexMisIDFV[pmtType] = new TH1D(Form("hvertexMisIDFV_pmtType%d",pmtType),"",20,0,200);
    hvertexMisIDFV_total[pmtType] = new TH1D(Form("hvertexMisIDFV_total_pmtType%d",pmtType),"",20,0,200);
    hvertexMisIDFV[pmtType]->GetXaxis()->SetTitle("FV distance to the wall");
    hvertexMisIDFV[pmtType]->GetYaxis()->SetTitle("Probability of migration from out to in FV");

    hvertexMisIDExtCut[pmtType] = new TH1D(Form("hvertexMisIDExtCut_pmtType%d",pmtType),"",20,0,2000);
    hvertexMisIDExtCut_total[pmtType] = new TH1D(Form("hvertexMisIDExtCut_total_pmtType%d",pmtType),"",20,0,2000);
    hvertexMisIDExtCut[pmtType]->GetXaxis()->SetTitle("ExtCut distance to the wall");
    hvertexMisIDExtCut[pmtType]->GetYaxis()->SetTitle("Probability of migration from out to in ExtCut");

    hvertex2D[pmtType] = new TH2D(Form("hvertex2D_pmtType%d",pmtType),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
    hvertex2D[pmtType]->GetXaxis()->SetTitle("x position (cm)");
    hvertex2D[pmtType]->GetYaxis()->SetTitle("y position (cm)");

    htruevertex2D[pmtType] = new TH2D(Form("htruevertex2D_pmtType%d",pmtType),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
    htruevertex2D[pmtType]->GetXaxis()->SetTitle("x position (cm)");
    htruevertex2D[pmtType]->GetYaxis()->SetTitle("y position (cm)");
  }

  bstree->SetBranchAddress("lf_time", &lf_time);
  bstree->SetBranchAddress("lf_vertex", lf_vertex);
  bstree->SetBranchAddress("lf_dir", lf_dir);
  
  bstree->SetBranchAddress("true_vertex_X", &true_vertex_X);
  bstree->SetBranchAddress("true_vertex_Y", &true_vertex_Y);
  bstree->SetBranchAddress("true_vertex_Z", &true_vertex_Z);
  bstree->SetBranchAddress("true_vertex_T", &true_vertex_T);

  bstree->SetBranchAddress("true_dir_X", &true_dir_X);
  bstree->SetBranchAddress("true_dir_Y", &true_dir_Y);
  bstree->SetBranchAddress("true_dir_Z", &true_dir_Z);
			     
  bstree->SetBranchAddress("lf_good", &lf_good);
  bstree->SetBranchAddress("lf_dirKS", &lf_dirKS);
  bstree->SetBranchAddress("ID_hits", &ID_hits);
  bstree->SetBranchAddress("lf_intime",&lf_intime);
  bstree->SetBranchAddress("ID_hits_200",&ID_hits_200);

  bstree->SetBranchAddress("eventId",&eventId);
  bstree->SetBranchAddress("triggerId",&triggerId);

  double mean_noise_BL=0;//34;
  double cut_BL=0;
  double mean_noise_mPMT=0;//15;
  double cut_mPMT=0;
  double hit_intime_BL=19;//11;//11;//2*17;
  double hit_intime_mPMT=2;//12;//12;
  bool FVCut = false;//true;//false;//true;//false;//true;
  bool ExtCut = false;//true;//false;//true;//false;//true;
  double FVvalue = 200;//in cm
  double ExtCutvalue = 1200;//in cm
  
  int nentries = bstree->GetEntries();
  cout<<"Number of entries = "<<nentries<<endl;
  int nhits_BL=0;//Number of hits only in BL PMTs
  double vertex_distance_BL=0;//Number of hits only in BL PMTs
  //double memory=0;     




  
  //######################################################################################################
  //######################################################################################################
  //############################## LOOP OVER EACH DAQ/TRUE VERTEX EVENT ################################## 
  //######################################################################################################
  //######################################################################################################


  for (int entry=0; entry<nentries; entry++) {
    bstree->GetEntry(entry);
    nEventsTotal[pmtT]++;

    //######################################################################################################
    //######################################################################################################
    //###################################### LOOP OVER EACH TRIGGERS #######################################
    //######################################################################################################
    //######################################################################################################

    //######################################################################################################
    // One DAQ/Sim event can contain several vertices, and several trigger per vertices. Here, we assume one
    // vertex per event, but allow several trigger for each of them.
    //######################################################################################################
    int ntriggers = true_vertex_X->size();
    //if(ntriggers!=1) cout<<"In event "<<eventId<<", there are "<<ntriggers<<" triggers"<<endl;
    for(int itr=0;itr<ntriggers;itr++){
      float mc_vertex[3] = {true_vertex_X->at(itr),true_vertex_Y->at(itr),true_vertex_Z->at(itr)};
      float mc_time = true_vertex_T->at(itr);
      float mc_dir[3] = {true_dir_X->at(itr),true_dir_Y->at(itr),true_dir_Z->at(itr)};

      //######################################################################################################
      //######################################################################################################
      //############################# Calculation of dWall and to wall #######################################
      //######################################################################################################
      //######################################################################################################
      
      double radius_true_vertex = TMath::Sqrt( pow(mc_vertex[0],2) + pow(mc_vertex[1],2) );
      //cout<<"Radius = "<<TMath::Sqrt( pow(mc_vertex[0],2) + pow(mc_vertex[1],2) )<<", and z position = "<<mc_vertex[2]<<endl;
      //dwall is the smallest of distance from vertex to top/bottom vs to edge
      double dwall = TMath::Min(tankHeight-TMath::Abs(mc_vertex[2]),TMath::Abs(radius_true_vertex-tankRadius));
      htruevertex2D[pmtT]->Fill(mc_vertex[0],mc_vertex[1]);
      hvertex2D[pmtT]->Fill(lf_vertex[0],lf_vertex[1]);
      double radius_rec_vertex = TMath::Sqrt( pow(lf_vertex[0],2) + pow(lf_vertex[1],2) );
      //dwall is the smallest of distance from vertex to top/bottom vs to edge
      double dwall_rec = TMath::Min(tankHeight-TMath::Abs(lf_vertex[2]),TMath::Abs(radius_rec_vertex-tankRadius));
      
      
      double towall=0;
      double stepSize=1;//in cm
      double position[3]={0.};
      double rad=0;double height=0;
      while(rad<tankRadius && TMath::Abs(height)<tankHeight){
	towall+=stepSize;
	for(int i=0;i<3;i++){
	  position[i]=mc_vertex[i]+mc_dir[i]*towall;
	}
	rad=TMath::Sqrt( pow(position[0],2) + pow(position[1],2) );
	height=position[2];
      }
      
      double pmt_to_vertex=0;
      rad=0;height=0;
      for(int i=0;i<3;i++) position[i]=0;
      while(rad<tankRadius && TMath::Abs(height)<tankHeight){
	pmt_to_vertex+=stepSize;
	for(int i=0;i<3;i++){
	  position[i]=mc_vertex[i]-mc_dir[i]*pmt_to_vertex;
	  //position[i]=mc_vertex[i]-mc_dir[i]*pmt_to_vertex;
	}
	rad=TMath::Sqrt( pow(position[0],2) + pow(position[1],2) );
	height=position[2];
      }
      
      double pmt_to_vertex_rec=0;
      rad=0;height=0;
      for(int i=0;i<3;i++) position[i]=0;
      while(rad<tankRadius && TMath::Abs(height)<tankHeight){
	pmt_to_vertex_rec+=stepSize;
	for(int i=0;i<3;i++){
	  position[i]=lf_vertex[i]-mc_dir[i]*pmt_to_vertex_rec;
	  //position[i]=mc_vertex[i]-mc_dir[i]*pmt_to_vertex;
	}
	rad=TMath::Sqrt( pow(position[0],2) + pow(position[1],2) );
	height=position[2];
      }
      
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      
      

      //######################################################################################################
      //######################################################################################################
      //#################################### Analysis cuts ###################################################
      //######################################################################################################
      //######################################################################################################

      if(FVCut && dwall<FVvalue) continue;
      if(ExtCut && pmt_to_vertex<ExtCutvalue) continue;
      nEvents[pmtT]++;

      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################

    

      //######################################################################################################
      //######################################################################################################
      //################################# Vertex resolution ##################################################
      //######################################################################################################
      //######################################################################################################

      //######################################################################################################
      // Simple vertex resolution
      //######################################################################################################
      double vertex_distance = sqrt(pow(lf_vertex[0]-mc_vertex[0],2) + pow(lf_vertex[1]-mc_vertex[1],2) + pow(lf_vertex[2]-mc_vertex[2],2));
      hvertex[pmtT]->Fill(vertex_distance);
      hradius[pmtT]->Fill(pow(radius_rec_vertex,2));
      hnhits[pmtT]->Fill(ID_hits);
      hvertexXnhits[pmtT]->Fill(vertex_distance,ID_hits);
      hdwallXnhits[pmtT]->Fill(dwall,ID_hits);
      htowallXnhits[pmtT]->Fill(towall,ID_hits);
      hvertexXnhits_intime[pmtT]->Fill(vertex_distance,lf_intime);
      nVertexFound[pmtT]++;

      //######################################################################################################
      // Vertex resolution along particle direction and transverse
      //######################################################################################################
      double vector_director[3];double vector_director_norm;
      for(int i=0;i<3;i++){
	vector_director[i]=lf_vertex[i]-mc_vertex[i];
	vector_director_norm+=pow(vector_director[i],2);
      }
      vector_director_norm = TMath::Sqrt(vector_director_norm);
      double vector_director_longitudinal = 0;//Component parralel to the direction of the particle
      double vector_director_ortho = 0;//Component ortho to the direction of the particle
      for(int i=0;i<3;i++){
	vector_director_longitudinal += vector_director[i]*mc_dir[i];
      }
      double angle=TMath::Abs(TMath::ACos(vector_director_longitudinal/vector_director_norm));
      //cout<<"angle = "<<angle*180/TMath::Pi()<<endl;
      vector_director_ortho = vector_director_norm*TMath::Sin(angle);
      //cout<<"Longitudinal vertex reso = "<<vector_director_longitudinal<<", time difference = "<<lf_vertex[3]-mc_vertex[3]<<endl;
      hvertexLongitudinal[pmtT]->Fill(TMath::Abs(vector_director_longitudinal));
      hvertexOrtho[pmtT]->Fill(vector_director_ortho);
      hvertexLongitudinalXTime[pmtT]->Fill(TMath::Abs(vector_director_longitudinal),lf_vertex[3]-mc_vertex[3]);

      //######################################################################################################
      // Vertex resolution vs dwall and towall
      //######################################################################################################
      hvertexXdwall[pmtT]->Fill(vertex_distance,dwall);
      hvertexXtowall[pmtT]->Fill(vertex_distance,towall);
      hdwallXdwallrec[pmtT]->Fill(dwall,(dwall_rec-dwall)/dwall);
     
      for(int ibinx=1;ibinx<=hvertexMisIDFV[pmtT]->GetNbinsX();ibinx++){
	double dwall_FV_definition=hvertexMisIDFV[pmtT]->GetBinLowEdge(ibinx)+hvertexMisIDFV[pmtT]->GetBinWidth(ibinx);
	if(dwall<dwall_FV_definition){//check only if the event is out of FV
	  bool misID=migrateFromOVtoFV(dwall, dwall_rec, FVvalue);
	  if(misID) hvertexMisIDFV[pmtT]->Fill(dwall_FV_definition-hvertexMisIDFV[pmtT]->GetBinWidth(ibinx));
	  hvertexMisIDFV_total[pmtT]->Fill(dwall_FV_definition-hvertexMisIDFV[pmtT]->GetBinWidth(ibinx));
	}
      }
    
      if(dwall<FVvalue){
	bool migration=migrateFromOVtoFV(dwall,dwall_rec,FVvalue);
	if(migration) nMisIDFV[pmtT]++;
	nMisIDFV_total[pmtT]++;
      }

      for(int ibinx=1;ibinx<=hvertexMisIDExtCut[pmtT]->GetNbinsX();ibinx++){
	double ExtCut_definition=hvertexMisIDExtCut[pmtT]->GetBinLowEdge(ibinx)+hvertexMisIDExtCut[pmtT]->GetBinWidth(ibinx);
	//cout<<dwall_ExtCut_definition<<endl;
	//if(pmt_to_vertex<ExtCut_definition){//check only if the event is out of ExtCut
	bool misID=migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_rec, ExtCut_definition);
	if(misID) hvertexMisIDExtCut[pmtT]->Fill(ExtCut_definition-hvertexMisIDExtCut[pmtT]->GetBinWidth(ibinx));
	hvertexMisIDExtCut_total[pmtT]->Fill(ExtCut_definition-hvertexMisIDExtCut[pmtT]->GetBinWidth(ibinx));
	//cout<<"misID="<<misID<<endl;
	//}
      }

      //if(pmt_to_vertex<ExtCutvalue){
      bool migration=migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_rec, ExtCutvalue);
      //bool migration=(pmt_to_vertex>ExtCutvalue);
      if(migration) nMisIDExtCut[pmtT]++;
      nMisIDExtCut_total[pmtT]++;

      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################

    

      //######################################################################################################
      //######################################################################################################
      //################################# Direction resolution ###############################################
      //######################################################################################################
      //######################################################################################################

      //######################################################################################################
      // Simple direction resolution
      //######################################################################################################
      double direction_angle = TMath::ACos(lf_dir[0]*mc_dir[0]+lf_dir[1]*mc_dir[1]+lf_dir[2]*mc_dir[2])*180/TMath::Pi();//Scalar product between true and rec dir.
      hdirection[pmtT]->Fill(direction_angle);
    
    
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################

    
    }

    //######################################################################################################
    //######################################################################################################
    //############################### END OF LOOP OVER EACH TRIGGERS #######################################
    //######################################################################################################
    //######################################################################################################

    
  }
  
  //######################################################################################################
  //######################################################################################################
  //########################## END OF LOOP OVER EACH DAQ/TRUE VERTEX EVENT ############################### 
  //######################################################################################################
  //######################################################################################################



  
 double misIDExtCut[nPMTtypes];
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    efficiency[pmtType] = (nEvents[pmtType]==0?0:nVertexFound[pmtType]/nEvents[pmtType]);
    misIDFV[pmtType] = (nEvents[pmtType]==0?0:nMisIDFV[pmtType]/nMisIDFV_total[pmtType]);
    //misIDFV[pmtType] = (nEvents[pmtType]==0?0:nMisIDExtCut[pmtType]/nMisIDExtCut_total[pmtType]);
    cout<<"Efficiency in loop : "<< efficiency[pmtType] << " and misID = "<<misIDFV[pmtType]<<endl;
    misIDExtCut[pmtType] = (nEvents[pmtType]==0?0:nMisIDExtCut[pmtType]/nMisIDExtCut_total[pmtType]);
    //misIDExtCut[pmtType] = (nEvents[pmtType]==0?0:nMisIDExtCut[pmtType]/nMisIDExtCut_total[pmtType]);
    cout<<"Efficiency in loop : "<< efficiency[pmtType] << " and misID = "<<misIDExtCut[pmtType]<<endl;
    cout<<"Number of events passing the FV and other cuts = "<<nEvents[pmtType]/nEventsTotal[pmtType]<<endl;
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    /*
    double integralSoFar = 0;
    double integralTotal = hvertex[pmtType]->Integral();
    double integral = 0;
    int ibin1D=1;
    while(integral<CredibleInterval && ibin1D<1e4){
      integralSoFar+=hvertex[pmtType]->GetBinContent(ibin1D);
      integral=integralSoFar/(integralTotal==0?1:integralTotal);
      ibin1D++;
      //cout<<integral<<endl;
    }
    vertex_resolution[pmtType] = hvertex[pmtType]->GetBinCenter(ibin1D);
*/
    vertex_resolution[pmtType] = vRes(hvertex[pmtType]);
    double error = vResAndError(hvertex[pmtType],CredibleInterval);
    //vResAndError
  }

  TCanvas * cvertex2D = new TCanvas("cvertex2D","");
  cvertex2D->Divide(2);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertex2D->cd(pmtType+1);
    hvertex2D[pmtType]->Draw("colz");
  }
  cvertex2D->SaveAs("plots/vertex2D.eps");
  

  TCanvas * ctruevertex2D = new TCanvas("ctruevertex2D","");
  ctruevertex2D->Divide(2);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    ctruevertex2D->cd(pmtType+1);
    htruevertex2D[pmtT]->Draw("colz");
  }
  ctruevertex2D->SaveAs("plots/truevertex2D.eps");

  
  TCanvas * cvertex = new TCanvas("cvertex","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hvertex[pmtType]->Rebin(10);
    hvertex[pmtType]->GetXaxis()->SetRangeUser(0,300);
    if(pmtType==0){
      hvertex[pmtType]->Scale(1./hvertex[pmtType]->Integral());
      hvertex[pmtType]->SetLineColor(kBlue);
      hvertex[pmtType]->Draw();
    }
    else{
      hvertex[pmtType]->Scale(1./hvertex[pmtType]->Integral());
      hvertex[pmtType]->SetLineColor(kRed);
      hvertex[pmtType]->Draw("same");
    }
    hvertex[pmtType]->SaveAs(Form("root_output/hvertex_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cvertex->SaveAs("plots/vertex.eps");

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexLongitudinal = new TCanvas("cvertexLongitudinal","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hvertexLongitudinal[pmtType]->Rebin(10);
    hvertexLongitudinal[pmtType]->GetXaxis()->SetRangeUser(0,300);
    if(pmtType==0){
      hvertexLongitudinal[pmtType]->SetLineColor(kBlue);
      hvertexLongitudinal[pmtType]->Draw();
    }
    else{
      hvertexLongitudinal[pmtType]->SetLineColor(kRed);
      hvertexLongitudinal[pmtType]->Draw("same");
    }
  }
  cvertexLongitudinal->SaveAs("plots/vertexLongitudinal.eps");
  
  TCanvas * cvertexOrtho = new TCanvas("cvertexOrtho","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hvertexOrtho[pmtType]->Rebin(10);
    hvertexOrtho[pmtType]->GetXaxis()->SetRangeUser(0,300);
    if(pmtType==0){
      hvertexOrtho[pmtType]->SetLineColor(kBlue);
      hvertexOrtho[pmtType]->Draw();
    }
    else{
      hvertexOrtho[pmtType]->SetLineColor(kRed);
      hvertexOrtho[pmtType]->Draw("same");
    }
  }
  cvertexOrtho->SaveAs("plots/vertexOrtho.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexLongitudinalXTime = new TCanvas("cvertexLongitudinalXTime","");
  cvertexLongitudinalXTime->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexLongitudinalXTime->cd(pmtType+1);
    hvertexLongitudinalXTime[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////



  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall = new TCanvas("cvertexXdwall","");
  cvertexXdwall->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexXdwall->cd(pmtType+1);
    //hvertexXdwall[pmtType]->Rebin2D(1,;
    hvertexXdwall[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall = new TCanvas("cvertexXtowall","");
  cvertexXtowall->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexXtowall->cd(pmtType+1);
    //hvertexXtowall[pmtType]->Rebin2D(1,;
    hvertexXtowall[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall1DAll = new TCanvas("cvertexXdwall1DAll","");
  const int nSlicesdwall = nbin_distance;//16;//hvertexXdwall[0]->GetNbinsX();
  cvertexXdwall1DAll->Divide(TMath::Sqrt(nSlicesdwall),TMath::Sqrt(nSlicesdwall));
  TH1D * hvertexXdwall1D[nPMTtypes][nSlicesdwall];
  double dwallPosition[nPMTtypes][nSlicesdwall];
  double resolution68[nPMTtypes][nSlicesdwall];
  double error_resolution68[nPMTtypes][nSlicesdwall];
  double error_x[nPMTtypes][nSlicesdwall];
  TGraphErrors * gdwall[nPMTtypes];
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    for(int ibinx=1;ibinx<=nSlicesdwall;ibinx++){
      hvertexXdwall1D[pmtType][ibinx-1] = (TH1D*) hvertexXdwall[pmtType]->ProjectionX(Form("hvertexXdwall1D%d_%d",pmtType,ibinx-1),ibinx,ibinx);

      dwallPosition[pmtType][ibinx-1]=hvertexXdwall[pmtType]->GetYaxis()->GetBinCenter(ibinx);
      resolution68[pmtType][ibinx-1]=vRes(hvertexXdwall1D[pmtType][ibinx-1]);
      error_resolution68[pmtType][ibinx-1]=vResAndError(hvertexXdwall1D[pmtType][ibinx-1],CredibleInterval);
      //cout<<"Error = "<<error_resolution68[pmtType][ibinx-1]<<endl;
      error_x[pmtType][ibinx-1]=0;
      //hvertexXdwall1D[pmtType][ibinx-1]->Scale(1./hvertexXdwall1D[pmtType][ibinx-1]->Integral());
      hvertexXdwall1D[pmtType][ibinx-1]->GetXaxis()->SetTitle(hvertexXdwall[pmtType]->GetYaxis()->GetTitle());
      if(pmtType==0){
	cvertexXdwall1DAll->cd(ibinx);
	hvertexXdwall1D[pmtType][ibinx-1]->SetLineColor(kBlue);
	hvertexXdwall1D[pmtType][ibinx-1]->Draw();
      }
      else{
	cvertexXdwall1DAll->cd(ibinx);
	hvertexXdwall1D[pmtType][ibinx-1]->SetLineColor(kRed);
	hvertexXdwall1D[pmtType][ibinx-1]->Draw("same");
      }
      /*
      double integralSoFar = 0;
      double integralTotal = hvertexXdwall1D[pmtType][ibinx-1]->Integral();
      double integral = 0;
      int ibin1D=1;
      while(integral<CredibleInterval && ibin1D<1e4){
	integralSoFar+=hvertexXdwall1D[pmtType][ibinx-1]->GetBinContent(ibin1D);
	integral=integralSoFar/(integralTotal==0?1:integralTotal);
	//cout<<"bin="<<ibin1D<<", integral="<<integral<<endl;
	ibin1D++;
	//cout<<integral<<endl;
      }
      dwallPosition[pmtType][ibinx-1]=hvertexXdwall[pmtType]->GetYaxis()->GetBinCenter(ibinx);
      resolution68[pmtType][ibinx-1]=hvertexXdwall1D[pmtType][ibinx-1]->GetBinCenter(ibin1D);
      */
      //cout<<"dwall="<<dwallPosition[pmtType][ibinx-1]<<", resolution="<<resolution68[pmtType][ibinx-1]<<endl;
    }
    gdwall[pmtType] = new TGraphErrors(nSlicesdwall,dwallPosition[pmtType],resolution68[pmtType],error_x[pmtType],error_resolution68[pmtType]);
    gdwall[pmtType]->GetXaxis()->SetTitle(hvertexXdwall[pmtType]->GetYaxis()->GetTitle());
    gdwall[pmtType]->GetYaxis()->SetTitle(hvertexXdwall[pmtType]->GetXaxis()->GetTitle());
    gdwall[pmtType]->SaveAs(Form("plots/gdwall_pmtType%d_%s.root",pmtType,file->GetName()));
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall1D = new TCanvas("cvertexXdwall1D","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    gdwall[pmtType]->SetMarkerStyle(20);
    gdwall[pmtType]->GetXaxis()->SetRangeUser(min_distance,max_distance);
    gdwall[pmtType]->GetYaxis()->SetRangeUser(50,100);
    if(pmtType==0){
      gdwall[pmtType]->SetLineColor(kBlue);
      gdwall[pmtType]->SetMarkerColor(kBlue);
      gdwall[pmtType]->Draw("ALP");
    }
    else{
      gdwall[pmtType]->SetMarkerColor(kRed);
      gdwall[pmtType]->SetLineColor(kRed);
      gdwall[pmtType]->Draw("PLsame");
    }
    gdwall[pmtType]->SaveAs(Form("plots/gdwall_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cvertexXdwall1D->SaveAs("plots/vertexXdwall.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////




  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall1DAll = new TCanvas("cvertexXtowall1DAll","");
  const int nSlicestowall = nbin_distance;//16;//hvertexXtowall[0]->GetNbinsX();
  cvertexXtowall1DAll->Divide(TMath::Sqrt(nSlicestowall),TMath::Sqrt(nSlicestowall));
  TH1D * hvertexXtowall1D[nPMTtypes][nSlicestowall];
  double towallPosition[nPMTtypes][nSlicestowall];
  double resolution68_towall[nPMTtypes][nSlicestowall];
  double error_resolution68_towall[nPMTtypes][nSlicestowall];
  double error_x_towall[nPMTtypes][nSlicestowall];
  TGraphErrors * gtowall[nPMTtypes];
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    for(int ibinx=1;ibinx<=nSlicestowall;ibinx++){
      hvertexXtowall1D[pmtType][ibinx-1] = (TH1D*) hvertexXtowall[pmtType]->ProjectionX(Form("hvertexXtowall1D%d_%d",pmtType,ibinx-1),ibinx,ibinx);

      towallPosition[pmtType][ibinx-1]=hvertexXtowall[pmtType]->GetYaxis()->GetBinCenter(ibinx);
      resolution68_towall[pmtType][ibinx-1]=vRes(hvertexXtowall1D[pmtType][ibinx-1]);
      error_resolution68_towall[pmtType][ibinx-1]=vResAndError(hvertexXtowall1D[pmtType][ibinx-1],CredibleInterval);
      //cout<<"Error = "<<error_resolution68_towall[pmtType][ibinx-1]<<endl;
      error_x_towall[pmtType][ibinx-1]=0;
      //hvertexXtowall1D[pmtType][ibinx-1]->Scale(1./hvertexXtowall1D[pmtType][ibinx-1]->Integral());
      hvertexXtowall1D[pmtType][ibinx-1]->GetXaxis()->SetTitle(hvertexXtowall[pmtType]->GetYaxis()->GetTitle());
      if(pmtType==0){
	cvertexXtowall1DAll->cd(ibinx);
	hvertexXtowall1D[pmtType][ibinx-1]->SetLineColor(kBlue);
	hvertexXtowall1D[pmtType][ibinx-1]->Draw();
      }
      else{
	cvertexXtowall1DAll->cd(ibinx);
	hvertexXtowall1D[pmtType][ibinx-1]->SetLineColor(kRed);
	hvertexXtowall1D[pmtType][ibinx-1]->Draw("same");
      }
      /*
      double integralSoFar = 0;
      double integralTotal = hvertexXtowall1D[pmtType][ibinx-1]->Integral();
      double integral = 0;
      int ibin1D=1;
      while(integral<CredibleInterval && ibin1D<1e4){
	integralSoFar+=hvertexXtowall1D[pmtType][ibinx-1]->GetBinContent(ibin1D);
	integral=integralSoFar/(integralTotal==0?1:integralTotal);
	//cout<<"bin="<<ibin1D<<", integral="<<integral<<endl;
	ibin1D++;
	//cout<<integral<<endl;
      }
      towallPosition[pmtType][ibinx-1]=hvertexXtowall[pmtType]->GetYaxis()->GetBinCenter(ibinx);
      resolution68_towall[pmtType][ibinx-1]=hvertexXtowall1D[pmtType][ibinx-1]->GetBinCenter(ibin1D);
      */
      //cout<<"towall="<<towallPosition[pmtType][ibinx-1]<<", resolution="<<resolution68_towall[pmtType][ibinx-1]<<endl;
    }
    gtowall[pmtType] = new TGraphErrors(nSlicestowall,towallPosition[pmtType],resolution68_towall[pmtType],error_x_towall[pmtType],error_resolution68_towall[pmtType]);
    gtowall[pmtType]->GetXaxis()->SetTitle(hvertexXtowall[pmtType]->GetYaxis()->GetTitle());
    gtowall[pmtType]->GetYaxis()->SetTitle(hvertexXtowall[pmtType]->GetXaxis()->GetTitle());
    gtowall[pmtType]->SaveAs(Form("plots/gtowall_pmtType%d_%s.root",pmtType,file->GetName()));
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall1D = new TCanvas("cvertexXtowall1D","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    gtowall[pmtType]->SetMarkerStyle(20);
    gtowall[pmtType]->GetXaxis()->SetRangeUser(min_distance,max_distance);
    gtowall[pmtType]->GetYaxis()->SetRangeUser(50,100);
    if(pmtType==0){
      gtowall[pmtType]->SetLineColor(kBlue);
      gtowall[pmtType]->SetMarkerColor(kBlue);
      gtowall[pmtType]->Draw("ALP");
    }
    else{
      gtowall[pmtType]->SetMarkerColor(kRed);
      gtowall[pmtType]->SetLineColor(kRed);
      gtowall[pmtType]->Draw("PLsame");
    }
    gtowall[pmtType]->SaveAs(Form("plots/gtowall_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cvertexXtowall1D->SaveAs("plots/vertexXtowall.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


    TCanvas * cradius = new TCanvas("cradius","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hradius[pmtType]->Rebin(10);
    //hradius[pmtType]->GetXaxis()->SetRangeUser(0,300);
    if(pmtType==0){
      hradius[pmtType]->Scale(1./hradius[pmtType]->Integral());
      hradius[pmtType]->SetLineColor(kBlue);
      hradius[pmtType]->Draw();
    }
    else{
      hradius[pmtType]->Scale(1./hradius[pmtType]->Integral());
      hradius[pmtType]->SetLineColor(kRed);
      hradius[pmtType]->Draw("same");
    }
    hradius[pmtType]->SaveAs(Form("root_output/hradius_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cradius->SaveAs("plots/radius.eps");


  
    ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdwallXdwallrec = new TCanvas("cdwallXdwallrec","");
  cdwallXdwallrec->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cdwallXdwallrec->cd(pmtType+1);
    //hdwallXdwallrec[pmtType]->Rebin2D(1,;
    hdwallXdwallrec[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cnhits = new TCanvas("cnhits","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    if(pmtType==0){
      hnhits[pmtType]->SetLineColor(kBlue);
      hnhits[pmtType]->Draw();
    }
    else{
      hnhits[pmtType]->SetLineColor(kRed);
      hnhits[pmtType]->Draw("same");
    }
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXnhits = new TCanvas("cvertexXnhits","");
  cvertexXnhits->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexXnhits->cd(pmtType+1);
    //hvertexXnhits[pmtType]->GetYaxis()->Get;
    hvertexXnhits[pmtType]->Rebin2D(1,2);
    hvertexXnhits[pmtType]->Draw("colz");
  }

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdwallXnhits = new TCanvas("cdwallXnhits","");
  //cdwallXnhits->Divide(2,1);
  for(int pmtType=0;pmtType<1;pmtType++){
    //cdwallXnhits->cd(pmtType+1);
    //hdwallXnhits[pmtType]->GetYaxis()->Get;
    hdwallXnhits[pmtType]->Rebin2D(1,2);
    hdwallXnhits[pmtType]->Draw("colz");
  }

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * ctowallXnhits = new TCanvas("ctowallXnhits","");
  //ctowallXnhits->Divide(2,1);
  for(int pmtType=0;pmtType<1;pmtType++){
    //ctowallXnhits->cd(pmtType+1);
    //htowallXnhits[pmtType]->GetYaxis()->Get;
    htowallXnhits[pmtType]->Rebin2D(1,2);
    htowallXnhits[pmtType]->Draw("colz");
  }

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXnhits_intime = new TCanvas("cvertexXnhits_intime","");
  cvertexXnhits_intime->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexXnhits_intime->cd(pmtType+1);
    hvertexXnhits_intime[pmtType]->Rebin2D(10,1);
    hvertexXnhits_intime[pmtType]->GetYaxis()->SetRangeUser(0,10);
    hvertexXnhits_intime[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////



  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXgoodness = new TCanvas("cvertexXgoodness","");
  cvertexXgoodness->Divide(2,1);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    cvertexXgoodness->cd(pmtType+1);
    //hvertexXgoodness[pmtType]->GetYaxis()->Get;
    hvertexXgoodness[pmtType]->Rebin2D(100,2);
    hvertexXgoodness[pmtType]->GetYaxis()->SetRangeUser(0.,1.);
    hvertexXgoodness[pmtType]->Draw("colz");
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


    ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    /*
    double integralSoFar = 0;
    double integralTotal = hdirection[pmtType]->Integral();
    double integral = 0;
    int ibin1D=1;
    while(integral<CredibleInterval && ibin1D<1e4){
      integralSoFar+=hdirection[pmtType]->GetBinContent(ibin1D);
      integral=integralSoFar/(integralTotal==0?1:integralTotal);
      ibin1D++;
      //cout<<integral<<endl;
    }
    cout<<"direction resolution ="<<hdirection[pmtType]->GetBinCenter(ibin1D)<<endl;*/
    //cout<<"direction resolution ="<<vRes(hdirection[pmtType])<<endl;    
    //direction_resolution[pmtType] = hdirection[pmtType]->GetBinCenter(ibin1D);
  }

  TCanvas * cdirection = new TCanvas("cdirection","");
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hdirection[pmtType]->Rebin(2);
    hdirection[pmtType]->GetXaxis()->SetRangeUser(0,100);
    if(pmtType==0){
      hdirection[pmtType]->Scale(1./hdirection[pmtType]->Integral());
      hdirection[pmtType]->SetLineColor(kBlue);
      hdirection[pmtType]->Draw();
    }
    else{
      hdirection[pmtType]->Scale(1./hdirection[pmtType]->Integral());
      hdirection[pmtType]->SetLineColor(kRed);
      hdirection[pmtType]->Draw("same");
    }
  }
  cdirection->SaveAs("plots/direction.eps");


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexMisIDFV = new TCanvas("cvertexMisIDFV","");
  gStyle->SetOptStat(0.);
  cvertexMisIDFV->SetGridx(true);
  cvertexMisIDFV->SetGridy(true);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hvertexMisIDFV[pmtType]->Divide(hvertexMisIDFV_total[pmtType]);
    hvertexMisIDFV[pmtType]->SetLineWidth(2);
    if(pmtType==0){
      hvertexMisIDFV[pmtType]->SetLineColor(kBlue);
      hvertexMisIDFV[pmtType]->GetXaxis()->SetRangeUser(0,200); 
      hvertexMisIDFV[pmtType]->GetYaxis()->SetRangeUser(0,1.); 
      hvertexMisIDFV[pmtType]->Draw();
    }
    else{
      hvertexMisIDFV[pmtType]->SetLineColor(kRed);
      hvertexMisIDFV[pmtType]->Draw("same");
    }
    hvertexMisIDFV[pmtType]->SaveAs(Form("root_output/hvertexMisIDFV_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cvertexMisIDFV->SaveAs("plots/vertexMisIDFV.eps");


  TCanvas * cvertexMisIDExtCut = new TCanvas("cvertexMisIDExtCut","");
  gStyle->SetOptStat(0.);
  cvertexMisIDExtCut->SetGridx(true);
  cvertexMisIDExtCut->SetGridy(true);
  for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
    hvertexMisIDExtCut[pmtType]->Divide(hvertexMisIDExtCut_total[pmtType]);
    hvertexMisIDExtCut[pmtType]->SetLineWidth(2);
    if(pmtType==0){
      hvertexMisIDExtCut[pmtType]->SetLineColor(kBlue);
      hvertexMisIDExtCut[pmtType]->GetXaxis()->SetRangeUser(0,2000); 
      hvertexMisIDExtCut[pmtType]->GetYaxis()->SetRangeUser(0,1.); 
      hvertexMisIDExtCut[pmtType]->Draw();
    }
    else{
      hvertexMisIDExtCut[pmtType]->SetLineColor(kRed);
      hvertexMisIDExtCut[pmtType]->Draw("same");
    }
    hvertexMisIDExtCut[pmtType]->SaveAs(Form("root_output/hvertexMisIDExtCut_pmttype%d_%s.root",pmtType,file->GetName()));
  }
  cvertexMisIDExtCut->SaveAs("plots/vertexMisIDExtCut.eps");

  cout<<"End of this file"<<endl;
  
  if(!monoEnergy){
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      delete hvertexXgoodness[pmtType];
      delete hvertexXnhits[pmtType];
      delete hvertexXdwall[pmtType];
      delete hvertexXtowall[pmtType];
      delete hvertexLongitudinalXTime[pmtType];delete hvertexLongitudinal[pmtType];delete hvertexOrtho[pmtType];
      delete hvertex[pmtType];delete hdirection[pmtType];delete hnhits[pmtType];
      delete hdwallXdwallrec[pmtType];
      delete gdwall[pmtType];
      delete gtowall[pmtType];
    }
    delete cvertexLongitudinal;
    delete cvertexOrtho;
    delete cvertex;
    delete cvertexXnhits;
    delete cnhits;
    delete cvertexXdwall1D;
    delete cvertexXdwall1DAll;
    delete cvertexXdwall;
    delete cvertexLongitudinalXTime;
  }
  return 0;
}

int main(int argc, char ** argv){
  //void bonsaiOutputAnalysisHybrid(){
  char filename[256];
  //sprintf(filename,"30PCmPMT0peCut");
  //sprintf(filename,"40PCBALOnly");
  //sprintf(filename,"50HzmPMT");
  //sprintf(filename,"25PCStandard");
  //sprintf(filename,"25PCStandard_2pecut");
  //sprintf(filename,"40PC_0pecut");
  //sprintf(filename,"30PC_2pecut");
  //sprintf(filename,"25PC_50Hz");
  //sprintf(filename,"MINUITbased_BLPDF");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_cutHits");
  //sprintf(filename,"Test");
  //sprintf(filename,"Test");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_40pc");
  //sprintf(filename,"MINUITbased_10PCmPMTPDF_faster_cutHits");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_6MeV");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_halfDR");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_10kmpmt_halfDR");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_40pcblonly_halfDR");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_halfDRonlyforBL");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_10kmpmt_halfDRonlyforBL");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_40pcblonly_halfMasked");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_10kmpmt_halfMasked");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_10kmpmt");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_0.8masked");
  //sprintf(filename,"MINUITbased_mPMTPDF_faster_directionality_0.4masked");
  sprintf(filename,"LEAF_5kmPMT");
  //sprintf(filename,"LEAF_40PCBL");
  //sprintf(filename,"LEAF_3kmPMT");
  //sprintf(filename,"test");
  //double PC = 40;PC=atof(argv[1]);
  //double DR = 4.2;DR=atof(argv[2]);
  double PC = 20;
  double DR = 1.0;//4.2;//8.4;
  //sprintf(filename,"LEAF_%.0fPC_%.1fkHz",PC,DR);
  //sprintf(filename,"LEAF_hybridPMT");
  gStyle->SetOptFit(1);
  TApplication theApp("theApp", &argc, argv);
  
  gStyle->SetOptTitle(0);
  const int nPMTtypes=1;
  double * vertex_resolution = new double[nPMTtypes];
  double * efficiency = new double[nPMTtypes];
  double * migrationProbability = new double[nPMTtypes];

  
  bool monoNRJ=false;//true;//false;//true;//false;
  bool energyOnly=false;//true;//true;
  
  if(monoNRJ){
    vertexResolution(3,vertex_resolution,efficiency,migrationProbability,nPMTtypes,monoNRJ,PC,DR);
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      cout<<"Efficiency ="<<efficiency[pmtType]<<endl;
      cout<<"vertex Resolution ="<<vertex_resolution[pmtType]<<endl;
    }
  }
  
  else if(energyOnly){    
    const int nNRJ = 7;double event_energy[nNRJ] = {3,4,5,6,8,10,15};
    double vertex_resolution_energy[nPMTtypes][nNRJ];
    double efficiency_energy[nPMTtypes][nNRJ];
    double migrationProbability_energy[nPMTtypes][nNRJ];
    
    for(int ie=0;ie<nNRJ;ie++){   
      vertexResolution(event_energy[ie],vertex_resolution,efficiency,migrationProbability,nPMTtypes,monoNRJ,PC,DR);
      for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	cout<<"PMT type ="<<pmtType<<", Efficiency ="<<efficiency[pmtType]<<", vertex Resolution ="<<vertex_resolution[pmtType]<<endl;
	vertex_resolution_energy[pmtType][ie] = vertex_resolution[pmtType];
	efficiency_energy[pmtType][ie] = efficiency[pmtType];
	migrationProbability_energy[pmtType][ie] = migrationProbability[pmtType];
      }
    }
    
    TCanvas * cEnergy = new TCanvas("cEnergy","");
    TGraph * genergy[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      genergy[pmtType] = new TGraph(nNRJ,event_energy,vertex_resolution_energy[pmtType]);
      genergy[pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
      genergy[pmtType]->GetYaxis()->SetTitle("Vertex resolution (cm)");
      genergy[pmtType]->GetXaxis()->SetRangeUser(3,15);
      genergy[pmtType]->GetYaxis()->SetRangeUser(0,300);
      genergy[pmtType]->SetMarkerStyle(20);
      if(pmtType==0){
	genergy[pmtType]->SetLineColor(kBlue);
	genergy[pmtType]->SetMarkerColor(kBlue);
	genergy[pmtType]->Draw("ALP");
      }
      else{
	genergy[pmtType]->SetLineColor(kRed);
	genergy[pmtType]->SetMarkerColor(kRed);
	genergy[pmtType]->Draw("sameLP");
      }
      genergy[pmtType]->SaveAs(Form("plots/genergy_pmtType%d_%s.root",pmtType,filename));
    }
    cEnergy->SaveAs("plots/vertexEnergy.eps");
    
    TCanvas * cEfficiency = new TCanvas("cEfficiency","");
    TGraph * gefficiency[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      gefficiency[pmtType] = new TGraph(nNRJ,event_energy,efficiency_energy[pmtType]);
      gefficiency[pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
      gefficiency[pmtType]->GetYaxis()->SetTitle("Signal efficiency");
      gefficiency[pmtType]->GetXaxis()->SetRangeUser(3,15);
      gefficiency[pmtType]->GetYaxis()->SetRangeUser(0,1.05);
      gefficiency[pmtType]->SetMarkerStyle(20);
      if(pmtType==0){
	gefficiency[pmtType]->SetLineColor(kBlue);
	gefficiency[pmtType]->SetMarkerColor(kBlue);
	gefficiency[pmtType]->Draw("ALP");
      }
      else{
	gefficiency[pmtType]->SetLineColor(kRed);
	gefficiency[pmtType]->SetMarkerColor(kRed);
	gefficiency[pmtType]->Draw("sameLP");
      }
      gefficiency[pmtType]->SaveAs(Form("plots/gefficiency_pmtType%d_%s.root",pmtType,filename));
    }
    cEfficiency->SaveAs("plots/vertexEfficiency.eps");    


    TCanvas * cMigrationProbability_mono = new TCanvas("cMigrationProbability_mono","");
    TGraph * gmigrationProbability[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      gmigrationProbability[pmtType] = new TGraph(nNRJ,event_energy,migrationProbability_energy[pmtType]);
      gmigrationProbability[pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
      gmigrationProbability[pmtType]->GetYaxis()->SetTitle("Signal migrationProbability");
      gmigrationProbability[pmtType]->GetXaxis()->SetRangeUser(3,15);
      gmigrationProbability[pmtType]->GetYaxis()->SetRangeUser(0,1.05);
      gmigrationProbability[pmtType]->SetMarkerStyle(20);
      if(pmtType==0){
	gmigrationProbability[pmtType]->SetLineColor(kBlue);
	gmigrationProbability[pmtType]->SetMarkerColor(kBlue);
	gmigrationProbability[pmtType]->Draw("ALP");
      }
      else{
	gmigrationProbability[pmtType]->SetLineColor(kRed);
	gmigrationProbability[pmtType]->SetMarkerColor(kRed);
	gmigrationProbability[pmtType]->Draw("sameLP");
      }
      gmigrationProbability[pmtType]->SaveAs(Form("plots/gmigrationProbability_pmtType%d_%s.root",pmtType,filename));
    }
    cMigrationProbability_mono->SaveAs("plots/vertexMigrationProbability.eps");    
  }


  else{
    const int nNRJ = 7;double event_energy[nNRJ] = {3,4,5,6,8,10,15};
    const int nPC = 3;double event_pc[nPC] = {20.,30.,40.};
    const int nDR = 4;double event_dr[nDR] = {4.2,6.0,8.4,10.};
    int color[nPC] = {600,416+2,1};//,632,416+2};
    int style[nDR] = {10,7,2,1};

    //const int nNRJ = 7;double event_energy[nNRJ] = {3,4,5,6,8,10,15};
    //const int nPC = 1;double event_pc[nPC] = {20.};//,30.,40.};
    //const int nDR = 4;double event_dr[nDR] = {4.2,6.0,8.4,10.};
    //int color[nPC] = {600};//,632,416+2};
    //int style[nDR] = {10,7,2,1};

    //const int nPC = 1;double event_pc[nPC] = {20.};
    //int color[nDR] = {1,600,416+2,632};//,632,416+2};
    //int style[nPC] = {10,2,1};
    //int color[nPC] = {2,600,632,416+2};
    //int color[nPC] = {2,600,632,416+2};
    //int color[nPC] = {2};
    double vertex_resolution_energy[nPC][nDR][nPMTtypes][nNRJ];
    double signal_to_noise[nPC][nDR][nPMTtypes][nNRJ];
    double efficiency_energy[nPC][nDR][nPMTtypes][nNRJ];
    double migrationProbability_energy[nPC][nDR][nPMTtypes][nNRJ];
    dataID.clear();
    int counter=0;
    bool minimize=false;
    

    if(minimize){
      /*      
      for(int ipc=0;ipc<nPC;ipc++){   
	for(int idr=0;idr<nDR;idr++){   
	  for(int ie=0;ie<nNRJ;ie++){
	    vertexResolution(event_energy[ie],vertex_resolution,efficiency,nPMTtypes,monoNRJ,event_pc[ipc],event_dr[idr]);
	    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	      //cout<<"PMT type ="<<pmtType<<", Efficiency ="<<efficiency[pmtType]<<", vertex Resolution ="<<vertex_resolution[pmtType]<<endl;
	      vertex_resolution_energy[ipc][idr][pmtType][ie] = vertex_resolution[pmtType];
	      efficiency_energy[ipc][idr][pmtType][ie] = efficiency[pmtType];
	      //signal_to_noise[ipc][idr][pmtType][ie]=SignalToNoise(event_energy[ie],event_pc[ipc],event_dr[idr],0,1,-0.01);
	      //signal_to_noise[ipc][idr][pmtType][ie]=SignalToNoise(event_energy[ie],event_pc[ipc],event_dr[idr],0,0.0076,-0.00021);
	      array <double,ntuple> data={event_pc[ipc],event_dr[ipc],event_energy[ie],pmtType,vertex_resolution[pmtType]};
	      dataID.push_back(data);
	      cout<<"array <double,ntuple> data"<<counter<<"={"<<event_pc[ipc]<<","<<event_dr[idr]<<","<<event_energy[ie]<<","<<pmtType<<","<<vertex_resolution[pmtType]<<"};dataID.push_back(data"<<counter<<");"<<endl;
	      counter++;
	    }
	  }
      }
      }*/
array <double,ntuple> data0={20,4.2,3,0,296.5};dataID.push_back(data0);
array <double,ntuple> data1={20,4.2,3,1,295.5};dataID.push_back(data1);
array <double,ntuple> data2={20,4.2,4,0,153.5};dataID.push_back(data2);
array <double,ntuple> data3={20,4.2,4,1,154.5};dataID.push_back(data3);
array <double,ntuple> data4={20,4.2,5,0,112.5};dataID.push_back(data4);
array <double,ntuple> data5={20,4.2,5,1,112.5};dataID.push_back(data5);
array <double,ntuple> data6={20,4.2,6,0,93.5};dataID.push_back(data6);
array <double,ntuple> data7={20,4.2,6,1,93.5};dataID.push_back(data7);
array <double,ntuple> data8={20,4.2,8,0,73.5};dataID.push_back(data8);
array <double,ntuple> data9={20,4.2,8,1,73.5};dataID.push_back(data9);
array <double,ntuple> data10={20,4.2,10,0,62.5};dataID.push_back(data10);
array <double,ntuple> data11={20,4.2,10,1,62.5};dataID.push_back(data11);
array <double,ntuple> data12={20,4.2,15,0,54.5};dataID.push_back(data12);
array <double,ntuple> data13={20,4.2,15,1,54.5};dataID.push_back(data13);
array <double,ntuple> data14={20,6,3,0,373.5};dataID.push_back(data14);
array <double,ntuple> data15={20,6,3,1,370.5};dataID.push_back(data15);
array <double,ntuple> data16={20,6,4,0,178.5};dataID.push_back(data16);
array <double,ntuple> data17={20,6,4,1,178.5};dataID.push_back(data17);
array <double,ntuple> data18={20,6,5,0,127.5};dataID.push_back(data18);
array <double,ntuple> data19={20,6,5,1,127.5};dataID.push_back(data19);
array <double,ntuple> data20={20,6,6,0,98.5};dataID.push_back(data20);
array <double,ntuple> data21={20,6,6,1,98.5};dataID.push_back(data21);
array <double,ntuple> data22={20,6,8,0,79.5};dataID.push_back(data22);
array <double,ntuple> data23={20,6,8,1,79.5};dataID.push_back(data23);
array <double,ntuple> data24={20,6,10,0,67.5};dataID.push_back(data24);
array <double,ntuple> data25={20,6,10,1,67.5};dataID.push_back(data25);
array <double,ntuple> data26={20,6,15,0,55.5};dataID.push_back(data26);
array <double,ntuple> data27={20,6,15,1,55.5};dataID.push_back(data27);
array <double,ntuple> data28={20,8.4,3,0,480.5};dataID.push_back(data28);
array <double,ntuple> data29={20,8.4,3,1,480.5};dataID.push_back(data29);
array <double,ntuple> data30={20,8.4,4,0,203.5};dataID.push_back(data30);
array <double,ntuple> data31={20,8.4,4,1,203.5};dataID.push_back(data31);
array <double,ntuple> data32={20,8.4,5,0,138.5};dataID.push_back(data32);
array <double,ntuple> data33={20,8.4,5,1,138.5};dataID.push_back(data33);
array <double,ntuple> data34={20,8.4,6,0,108.5};dataID.push_back(data34);
array <double,ntuple> data35={20,8.4,6,1,108.5};dataID.push_back(data35);
array <double,ntuple> data36={20,8.4,8,0,82.5};dataID.push_back(data36);
array <double,ntuple> data37={20,8.4,8,1,82.5};dataID.push_back(data37);
array <double,ntuple> data38={20,8.4,10,0,70.5};dataID.push_back(data38);
array <double,ntuple> data39={20,8.4,10,1,70.5};dataID.push_back(data39);
array <double,ntuple> data40={20,8.4,15,0,57.5};dataID.push_back(data40);
array <double,ntuple> data41={20,8.4,15,1,57.5};dataID.push_back(data41);
array <double,ntuple> data42={20,10,3,0,582.5};dataID.push_back(data42);
array <double,ntuple> data43={20,10,3,1,587.5};dataID.push_back(data43);
array <double,ntuple> data44={20,10,4,0,232.5};dataID.push_back(data44);
array <double,ntuple> data45={20,10,4,1,233.5};dataID.push_back(data45);
array <double,ntuple> data46={20,10,5,0,149.5};dataID.push_back(data46);
array <double,ntuple> data47={20,10,5,1,149.5};dataID.push_back(data47);
array <double,ntuple> data48={20,10,6,0,115.5};dataID.push_back(data48);
array <double,ntuple> data49={20,10,6,1,115.5};dataID.push_back(data49);
array <double,ntuple> data50={20,10,8,0,86.5};dataID.push_back(data50);
array <double,ntuple> data51={20,10,8,1,86.5};dataID.push_back(data51);
array <double,ntuple> data52={20,10,10,0,71.5};dataID.push_back(data52);
array <double,ntuple> data53={20,10,10,1,71.5};dataID.push_back(data53);
array <double,ntuple> data54={20,10,15,0,59.5};dataID.push_back(data54);
array <double,ntuple> data55={20,10,15,1,59.5};dataID.push_back(data55);
array <double,ntuple> data56={30,4.2,3,0,171.5};dataID.push_back(data56);
array <double,ntuple> data57={30,4.2,3,1,171.5};dataID.push_back(data57);
array <double,ntuple> data58={30,4.2,4,0,109.5};dataID.push_back(data58);
array <double,ntuple> data59={30,4.2,4,1,109.5};dataID.push_back(data59);
array <double,ntuple> data60={30,4.2,5,0,82.5};dataID.push_back(data60);
array <double,ntuple> data61={30,4.2,5,1,82.5};dataID.push_back(data61);
array <double,ntuple> data62={30,4.2,6,0,72.5};dataID.push_back(data62);
array <double,ntuple> data63={30,4.2,6,1,72.5};dataID.push_back(data63);
array <double,ntuple> data64={30,4.2,8,0,59.5};dataID.push_back(data64);
array <double,ntuple> data65={30,4.2,8,1,59.5};dataID.push_back(data65);
array <double,ntuple> data66={30,4.2,10,0,52.5};dataID.push_back(data66);
array <double,ntuple> data67={30,4.2,10,1,52.5};dataID.push_back(data67);
array <double,ntuple> data68={30,4.2,15,0,46.5};dataID.push_back(data68);
array <double,ntuple> data69={30,4.2,15,1,46.5};dataID.push_back(data69);
array <double,ntuple> data70={30,6,3,0,207.5};dataID.push_back(data70);
array <double,ntuple> data71={30,6,3,1,206.5};dataID.push_back(data71);
array <double,ntuple> data72={30,6,4,0,119.5};dataID.push_back(data72);
array <double,ntuple> data73={30,6,4,1,119.5};dataID.push_back(data73);
array <double,ntuple> data74={30,6,5,0,93.5};dataID.push_back(data74);
array <double,ntuple> data75={30,6,5,1,93.5};dataID.push_back(data75);
array <double,ntuple> data76={30,6,6,0,75.5};dataID.push_back(data76);
array <double,ntuple> data77={30,6,6,1,75.5};dataID.push_back(data77);
array <double,ntuple> data78={30,6,8,0,63.5};dataID.push_back(data78);
array <double,ntuple> data79={30,6,8,1,63.5};dataID.push_back(data79);
array <double,ntuple> data80={30,6,10,0,55.5};dataID.push_back(data80);
array <double,ntuple> data81={30,6,10,1,55.5};dataID.push_back(data81);
array <double,ntuple> data82={30,6,15,0,47.5};dataID.push_back(data82);
array <double,ntuple> data83={30,6,15,1,47.5};dataID.push_back(data83);
array <double,ntuple> data84={30,8.4,3,0,252.5};dataID.push_back(data84);
array <double,ntuple> data85={30,8.4,3,1,251.5};dataID.push_back(data85);
array <double,ntuple> data86={30,8.4,4,0,135.5};dataID.push_back(data86);
array <double,ntuple> data87={30,8.4,4,1,135.5};dataID.push_back(data87);
array <double,ntuple> data88={30,8.4,5,0,97.5};dataID.push_back(data88);
array <double,ntuple> data89={30,8.4,5,1,97.5};dataID.push_back(data89);
array <double,ntuple> data90={30,8.4,6,0,82.5};dataID.push_back(data90);
array <double,ntuple> data91={30,8.4,6,1,82.5};dataID.push_back(data91);
array <double,ntuple> data92={30,8.4,8,0,67.5};dataID.push_back(data92);
array <double,ntuple> data93={30,8.4,8,1,66.5};dataID.push_back(data93);
array <double,ntuple> data94={30,8.4,10,0,58.5};dataID.push_back(data94);
array <double,ntuple> data95={30,8.4,10,1,58.5};dataID.push_back(data95);
array <double,ntuple> data96={30,8.4,15,0,48.5};dataID.push_back(data96);
array <double,ntuple> data97={30,8.4,15,1,48.5};dataID.push_back(data97);
array <double,ntuple> data98={30,10,3,0,313.5};dataID.push_back(data98);
array <double,ntuple> data99={30,10,3,1,314.5};dataID.push_back(data99);
array <double,ntuple> data100={30,10,4,0,152.5};dataID.push_back(data100);
array <double,ntuple> data101={30,10,4,1,152.5};dataID.push_back(data101);
array <double,ntuple> data102={30,10,5,0,107.5};dataID.push_back(data102);
array <double,ntuple> data103={30,10,5,1,107.5};dataID.push_back(data103);
array <double,ntuple> data104={30,10,6,0,87.5};dataID.push_back(data104);
array <double,ntuple> data105={30,10,6,1,87.5};dataID.push_back(data105);
array <double,ntuple> data106={30,10,8,0,68.5};dataID.push_back(data106);
array <double,ntuple> data107={30,10,8,1,68.5};dataID.push_back(data107);
array <double,ntuple> data108={30,10,10,0,58.5};dataID.push_back(data108);
array <double,ntuple> data109={30,10,10,1,58.5};dataID.push_back(data109);
array <double,ntuple> data110={30,10,15,0,50.5};dataID.push_back(data110);
array <double,ntuple> data111={30,10,15,1,50.5};dataID.push_back(data111);
array <double,ntuple> data112={40,4.2,3,0,131.5};dataID.push_back(data112);
array <double,ntuple> data113={40,4.2,3,1,130.5};dataID.push_back(data113);
array <double,ntuple> data114={40,4.2,4,0,87.5};dataID.push_back(data114);
array <double,ntuple> data115={40,4.2,4,1,87.5};dataID.push_back(data115);
array <double,ntuple> data116={40,4.2,5,0,69.5};dataID.push_back(data116);
array <double,ntuple> data117={40,4.2,5,1,69.5};dataID.push_back(data117);
array <double,ntuple> data118={40,4.2,6,0,61.5};dataID.push_back(data118);
array <double,ntuple> data119={40,4.2,6,1,61.5};dataID.push_back(data119);
array <double,ntuple> data120={40,4.2,8,0,51.5};dataID.push_back(data120);
array <double,ntuple> data121={40,4.2,8,1,51.5};dataID.push_back(data121);
array <double,ntuple> data122={40,4.2,10,0,46.5};dataID.push_back(data122);
array <double,ntuple> data123={40,4.2,10,1,46.5};dataID.push_back(data123);
array <double,ntuple> data124={40,4.2,15,0,42.5};dataID.push_back(data124);
array <double,ntuple> data125={40,4.2,15,1,42.5};dataID.push_back(data125);
array <double,ntuple> data126={40,6,3,0,156.5};dataID.push_back(data126);
array <double,ntuple> data127={40,6,3,1,156.5};dataID.push_back(data127);
array <double,ntuple> data128={40,6,4,0,96.5};dataID.push_back(data128);
array <double,ntuple> data129={40,6,4,1,96.5};dataID.push_back(data129);
array <double,ntuple> data130={40,6,5,0,77.5};dataID.push_back(data130);
array <double,ntuple> data131={40,6,5,1,77.5};dataID.push_back(data131);
array <double,ntuple> data132={40,6,6,0,63.5};dataID.push_back(data132);
array <double,ntuple> data133={40,6,6,1,63.5};dataID.push_back(data133);
array <double,ntuple> data134={40,6,8,0,54.5};dataID.push_back(data134);
array <double,ntuple> data135={40,6,8,1,54.5};dataID.push_back(data135);
array <double,ntuple> data136={40,6,10,0,48.5};dataID.push_back(data136);
array <double,ntuple> data137={40,6,10,1,48.5};dataID.push_back(data137);
array <double,ntuple> data138={40,6,15,0,43.5};dataID.push_back(data138);
array <double,ntuple> data139={40,6,15,1,43.5};dataID.push_back(data139);
array <double,ntuple> data140={40,8.4,3,0,184.5};dataID.push_back(data140);
array <double,ntuple> data141={40,8.4,3,1,183.5};dataID.push_back(data141);
array <double,ntuple> data142={40,8.4,4,0,107.5};dataID.push_back(data142);
array <double,ntuple> data143={40,8.4,4,1,107.5};dataID.push_back(data143);
array <double,ntuple> data144={40,8.4,5,0,81.5};dataID.push_back(data144);
array <double,ntuple> data145={40,8.4,5,1,81.5};dataID.push_back(data145);
array <double,ntuple> data146={40,8.4,6,0,69.5};dataID.push_back(data146);
array <double,ntuple> data147={40,8.4,6,1,69.5};dataID.push_back(data147);
array <double,ntuple> data148={40,8.4,8,0,58.5};dataID.push_back(data148);
array <double,ntuple> data149={40,8.4,8,1,58.5};dataID.push_back(data149);
array <double,ntuple> data150={40,8.4,10,0,51.5};dataID.push_back(data150);
array <double,ntuple> data151={40,8.4,10,1,51.5};dataID.push_back(data151);
array <double,ntuple> data152={40,8.4,15,0,44.5};dataID.push_back(data152);
array <double,ntuple> data153={40,8.4,15,1,44.5};dataID.push_back(data153);
array <double,ntuple> data154={40,10,3,0,219.5};dataID.push_back(data154);
array <double,ntuple> data155={40,10,3,1,219.5};dataID.push_back(data155);
array <double,ntuple> data156={40,10,4,0,118.5};dataID.push_back(data156);
array <double,ntuple> data157={40,10,4,1,118.5};dataID.push_back(data157);
array <double,ntuple> data158={40,10,5,0,86.5};dataID.push_back(data158);
array <double,ntuple> data159={40,10,5,1,86.5};dataID.push_back(data159);
array <double,ntuple> data160={40,10,6,0,73.5};dataID.push_back(data160);
array <double,ntuple> data161={40,10,6,1,73.5};dataID.push_back(data161);
array <double,ntuple> data162={40,10,8,0,59.5};dataID.push_back(data162);
array <double,ntuple> data163={40,10,8,1,59.5};dataID.push_back(data163);
array <double,ntuple> data164={40,10,10,0,52.5};dataID.push_back(data164);
array <double,ntuple> data165={40,10,10,1,52.5};dataID.push_back(data165);
array <double,ntuple> data166={40,10,15,0,45.5};dataID.push_back(data166);
array <double,ntuple> data167={40,10,15,1,45.5};dataID.push_back(data167);
      
    /*
      array <double,ntuple> data0={20,4.2,3,0,296.5};dataID.push_back(data0);
      array <double,ntuple> data1={20,4.2,3,1,295.5};dataID.push_back(data1);
      array <double,ntuple> data2={20,4.2,4,0,153.5};dataID.push_back(data2);
      array <double,ntuple> data3={20,4.2,4,1,154.5};dataID.push_back(data3);
      array <double,ntuple> data4={20,4.2,5,0,112.5};dataID.push_back(data4);
      array <double,ntuple> data5={20,4.2,5,1,112.5};dataID.push_back(data5);
      array <double,ntuple> data6={20,4.2,6,0,93.5};dataID.push_back(data6);
      array <double,ntuple> data7={20,4.2,6,1,93.5};dataID.push_back(data7);
      array <double,ntuple> data8={20,4.2,8,0,73.5};dataID.push_back(data8);
      array <double,ntuple> data9={20,4.2,8,1,73.5};dataID.push_back(data9);
      array <double,ntuple> data10={20,4.2,10,0,62.5};dataID.push_back(data10);
      array <double,ntuple> data11={20,4.2,10,1,62.5};dataID.push_back(data11);
      array <double,ntuple> data12={20,4.2,15,0,54.5};dataID.push_back(data12);
      array <double,ntuple> data13={20,4.2,15,1,54.5};dataID.push_back(data13);
      array <double,ntuple> data14={20,4.2,3,0,373.5};dataID.push_back(data14);
      array <double,ntuple> data15={20,4.2,3,1,370.5};dataID.push_back(data15);
      array <double,ntuple> data16={20,4.2,4,0,178.5};dataID.push_back(data16);
      array <double,ntuple> data17={20,4.2,4,1,178.5};dataID.push_back(data17);
      array <double,ntuple> data18={20,4.2,5,0,127.5};dataID.push_back(data18);
      array <double,ntuple> data19={20,4.2,5,1,127.5};dataID.push_back(data19);
      array <double,ntuple> data20={20,4.2,6,0,98.5};dataID.push_back(data20);
      array <double,ntuple> data21={20,4.2,6,1,98.5};dataID.push_back(data21);
      array <double,ntuple> data22={20,4.2,8,0,79.5};dataID.push_back(data22);
      array <double,ntuple> data23={20,4.2,8,1,79.5};dataID.push_back(data23);
      array <double,ntuple> data24={20,4.2,10,0,67.5};dataID.push_back(data24);
      array <double,ntuple> data25={20,4.2,10,1,67.5};dataID.push_back(data25);
      array <double,ntuple> data26={20,4.2,15,0,55.5};dataID.push_back(data26);
      array <double,ntuple> data27={20,4.2,15,1,55.5};dataID.push_back(data27);
      array <double,ntuple> data28={20,4.2,3,0,480.5};dataID.push_back(data28);
      array <double,ntuple> data29={20,4.2,3,1,480.5};dataID.push_back(data29);
      array <double,ntuple> data30={20,4.2,4,0,203.5};dataID.push_back(data30);
      array <double,ntuple> data31={20,4.2,4,1,203.5};dataID.push_back(data31);
      array <double,ntuple> data32={20,4.2,5,0,138.5};dataID.push_back(data32);
      array <double,ntuple> data33={20,4.2,5,1,138.5};dataID.push_back(data33);
      array <double,ntuple> data34={20,4.2,6,0,108.5};dataID.push_back(data34);
      array <double,ntuple> data35={20,4.2,6,1,108.5};dataID.push_back(data35);
      array <double,ntuple> data36={20,4.2,8,0,82.5};dataID.push_back(data36);
      array <double,ntuple> data37={20,4.2,8,1,82.5};dataID.push_back(data37);
      array <double,ntuple> data38={20,4.2,10,0,70.5};dataID.push_back(data38);
      array <double,ntuple> data39={20,4.2,10,1,70.5};dataID.push_back(data39);
      array <double,ntuple> data40={20,4.2,15,0,57.5};dataID.push_back(data40);
      array <double,ntuple> data41={20,4.2,15,1,57.5};dataID.push_back(data41);
      array <double,ntuple> data42={20,4.2,3,0,582.5};dataID.push_back(data42);
      array <double,ntuple> data43={20,4.2,3,1,587.5};dataID.push_back(data43);
      array <double,ntuple> data44={20,4.2,4,0,232.5};dataID.push_back(data44);
      array <double,ntuple> data45={20,4.2,4,1,233.5};dataID.push_back(data45);
      array <double,ntuple> data46={20,4.2,5,0,149.5};dataID.push_back(data46);
      array <double,ntuple> data47={20,4.2,5,1,149.5};dataID.push_back(data47);
      array <double,ntuple> data48={20,4.2,6,0,115.5};dataID.push_back(data48);
      array <double,ntuple> data49={20,4.2,6,1,115.5};dataID.push_back(data49);
      array <double,ntuple> data50={20,4.2,8,0,86.5};dataID.push_back(data50);
      array <double,ntuple> data51={20,4.2,8,1,86.5};dataID.push_back(data51);
      array <double,ntuple> data52={20,4.2,10,0,71.5};dataID.push_back(data52);
      array <double,ntuple> data53={20,4.2,10,1,71.5};dataID.push_back(data53);
      array <double,ntuple> data54={20,4.2,15,0,59.5};dataID.push_back(data54);
      array <double,ntuple> data55={20,4.2,15,1,59.5};dataID.push_back(data55);
      array <double,ntuple> data56={30,6,3,0,171.5};dataID.push_back(data56);
      array <double,ntuple> data57={30,6,3,1,171.5};dataID.push_back(data57);
      array <double,ntuple> data58={30,6,4,0,109.5};dataID.push_back(data58);
      array <double,ntuple> data59={30,6,4,1,109.5};dataID.push_back(data59);
      array <double,ntuple> data60={30,6,5,0,82.5};dataID.push_back(data60);
      array <double,ntuple> data61={30,6,5,1,82.5};dataID.push_back(data61);
      array <double,ntuple> data62={30,6,6,0,72.5};dataID.push_back(data62);
      array <double,ntuple> data63={30,6,6,1,72.5};dataID.push_back(data63);
      array <double,ntuple> data64={30,6,8,0,59.5};dataID.push_back(data64);
      array <double,ntuple> data65={30,6,8,1,59.5};dataID.push_back(data65);
      array <double,ntuple> data66={30,6,10,0,52.5};dataID.push_back(data66);
      array <double,ntuple> data67={30,6,10,1,52.5};dataID.push_back(data67);
      array <double,ntuple> data68={30,6,15,0,46.5};dataID.push_back(data68);
      array <double,ntuple> data69={30,6,15,1,46.5};dataID.push_back(data69);
      array <double,ntuple> data70={30,6,3,0,207.5};dataID.push_back(data70);
      array <double,ntuple> data71={30,6,3,1,206.5};dataID.push_back(data71);
      array <double,ntuple> data72={30,6,4,0,119.5};dataID.push_back(data72);
      array <double,ntuple> data73={30,6,4,1,119.5};dataID.push_back(data73);
      array <double,ntuple> data74={30,6,5,0,93.5};dataID.push_back(data74);
      array <double,ntuple> data75={30,6,5,1,93.5};dataID.push_back(data75);
      array <double,ntuple> data76={30,6,6,0,75.5};dataID.push_back(data76);
      array <double,ntuple> data77={30,6,6,1,75.5};dataID.push_back(data77);
      array <double,ntuple> data78={30,6,8,0,63.5};dataID.push_back(data78);
      array <double,ntuple> data79={30,6,8,1,63.5};dataID.push_back(data79);
      array <double,ntuple> data80={30,6,10,0,55.5};dataID.push_back(data80);
      array <double,ntuple> data81={30,6,10,1,55.5};dataID.push_back(data81);
      array <double,ntuple> data82={30,6,15,0,47.5};dataID.push_back(data82);
      array <double,ntuple> data83={30,6,15,1,47.5};dataID.push_back(data83);
      array <double,ntuple> data84={30,6,3,0,252.5};dataID.push_back(data84);
      array <double,ntuple> data85={30,6,3,1,251.5};dataID.push_back(data85);
      array <double,ntuple> data86={30,6,4,0,135.5};dataID.push_back(data86);
      array <double,ntuple> data87={30,6,4,1,135.5};dataID.push_back(data87);
      array <double,ntuple> data88={30,6,5,0,97.5};dataID.push_back(data88);
      array <double,ntuple> data89={30,6,5,1,97.5};dataID.push_back(data89);
      array <double,ntuple> data90={30,6,6,0,82.5};dataID.push_back(data90);
      array <double,ntuple> data91={30,6,6,1,82.5};dataID.push_back(data91);
      array <double,ntuple> data92={30,6,8,0,67.5};dataID.push_back(data92);
      array <double,ntuple> data93={30,6,8,1,66.5};dataID.push_back(data93);
      array <double,ntuple> data94={30,6,10,0,58.5};dataID.push_back(data94);
      array <double,ntuple> data95={30,6,10,1,58.5};dataID.push_back(data95);
      array <double,ntuple> data96={30,6,15,0,48.5};dataID.push_back(data96);
      array <double,ntuple> data97={30,6,15,1,48.5};dataID.push_back(data97);
      array <double,ntuple> data98={30,6,3,0,313.5};dataID.push_back(data98);
      array <double,ntuple> data99={30,6,3,1,314.5};dataID.push_back(data99);
      array <double,ntuple> data100={30,6,4,0,152.5};dataID.push_back(data100);
      array <double,ntuple> data101={30,6,4,1,152.5};dataID.push_back(data101);
      array <double,ntuple> data102={30,6,5,0,107.5};dataID.push_back(data102);
      array <double,ntuple> data103={30,6,5,1,107.5};dataID.push_back(data103);
      array <double,ntuple> data104={30,6,6,0,87.5};dataID.push_back(data104);
      array <double,ntuple> data105={30,6,6,1,87.5};dataID.push_back(data105);
      array <double,ntuple> data106={30,6,8,0,68.5};dataID.push_back(data106);
      array <double,ntuple> data107={30,6,8,1,68.5};dataID.push_back(data107);
      array <double,ntuple> data108={30,6,10,0,58.5};dataID.push_back(data108);
      array <double,ntuple> data109={30,6,10,1,58.5};dataID.push_back(data109);
      array <double,ntuple> data110={30,6,15,0,50.5};dataID.push_back(data110);
      array <double,ntuple> data111={30,6,15,1,50.5};dataID.push_back(data111);
      array <double,ntuple> data112={40,8.4,3,0,131.5};dataID.push_back(data112);
      array <double,ntuple> data113={40,8.4,3,1,130.5};dataID.push_back(data113);
      array <double,ntuple> data114={40,8.4,4,0,87.5};dataID.push_back(data114);
      array <double,ntuple> data115={40,8.4,4,1,87.5};dataID.push_back(data115);
      array <double,ntuple> data116={40,8.4,5,0,69.5};dataID.push_back(data116);
      array <double,ntuple> data117={40,8.4,5,1,69.5};dataID.push_back(data117);
      array <double,ntuple> data118={40,8.4,6,0,61.5};dataID.push_back(data118);
      array <double,ntuple> data119={40,8.4,6,1,61.5};dataID.push_back(data119);
      array <double,ntuple> data120={40,8.4,8,0,51.5};dataID.push_back(data120);
      array <double,ntuple> data121={40,8.4,8,1,51.5};dataID.push_back(data121);
      array <double,ntuple> data122={40,8.4,10,0,46.5};dataID.push_back(data122);
      array <double,ntuple> data123={40,8.4,10,1,46.5};dataID.push_back(data123);
      array <double,ntuple> data124={40,8.4,15,0,42.5};dataID.push_back(data124);
      array <double,ntuple> data125={40,8.4,15,1,42.5};dataID.push_back(data125);
      array <double,ntuple> data126={40,8.4,3,0,156.5};dataID.push_back(data126);
      array <double,ntuple> data127={40,8.4,3,1,156.5};dataID.push_back(data127);
      array <double,ntuple> data128={40,8.4,4,0,96.5};dataID.push_back(data128);
      array <double,ntuple> data129={40,8.4,4,1,96.5};dataID.push_back(data129);
      array <double,ntuple> data130={40,8.4,5,0,77.5};dataID.push_back(data130);
      array <double,ntuple> data131={40,8.4,5,1,77.5};dataID.push_back(data131);
      array <double,ntuple> data132={40,8.4,6,0,63.5};dataID.push_back(data132);
      array <double,ntuple> data133={40,8.4,6,1,63.5};dataID.push_back(data133);
      array <double,ntuple> data134={40,8.4,8,0,54.5};dataID.push_back(data134);
      array <double,ntuple> data135={40,8.4,8,1,54.5};dataID.push_back(data135);
      array <double,ntuple> data136={40,8.4,10,0,48.5};dataID.push_back(data136);
      array <double,ntuple> data137={40,8.4,10,1,48.5};dataID.push_back(data137);
      array <double,ntuple> data138={40,8.4,15,0,43.5};dataID.push_back(data138);
      array <double,ntuple> data139={40,8.4,15,1,43.5};dataID.push_back(data139);
      array <double,ntuple> data140={40,8.4,3,0,184.5};dataID.push_back(data140);
      array <double,ntuple> data141={40,8.4,3,1,183.5};dataID.push_back(data141);
      array <double,ntuple> data142={40,8.4,4,0,107.5};dataID.push_back(data142);
      array <double,ntuple> data143={40,8.4,4,1,107.5};dataID.push_back(data143);
      array <double,ntuple> data144={40,8.4,5,0,81.5};dataID.push_back(data144);
      array <double,ntuple> data145={40,8.4,5,1,81.5};dataID.push_back(data145);
      array <double,ntuple> data146={40,8.4,6,0,69.5};dataID.push_back(data146);
      array <double,ntuple> data147={40,8.4,6,1,69.5};dataID.push_back(data147);
      array <double,ntuple> data148={40,8.4,8,0,58.5};dataID.push_back(data148);
      array <double,ntuple> data149={40,8.4,8,1,58.5};dataID.push_back(data149);
      array <double,ntuple> data150={40,8.4,10,0,51.5};dataID.push_back(data150);
      array <double,ntuple> data151={40,8.4,10,1,51.5};dataID.push_back(data151);
      array <double,ntuple> data152={40,8.4,15,0,44.5};dataID.push_back(data152);
      array <double,ntuple> data153={40,8.4,15,1,44.5};dataID.push_back(data153);
      array <double,ntuple> data154={40,8.4,3,0,219.5};dataID.push_back(data154);
      array <double,ntuple> data155={40,8.4,3,1,219.5};dataID.push_back(data155);
      array <double,ntuple> data156={40,8.4,4,0,118.5};dataID.push_back(data156);
      array <double,ntuple> data157={40,8.4,4,1,118.5};dataID.push_back(data157);
      array <double,ntuple> data158={40,8.4,5,0,86.5};dataID.push_back(data158);
      array <double,ntuple> data159={40,8.4,5,1,86.5};dataID.push_back(data159);
      array <double,ntuple> data160={40,8.4,6,0,73.5};dataID.push_back(data160);
      array <double,ntuple> data161={40,8.4,6,1,73.5};dataID.push_back(data161);
      array <double,ntuple> data162={40,8.4,8,0,59.5};dataID.push_back(data162);
      array <double,ntuple> data163={40,8.4,8,1,59.5};dataID.push_back(data163);
      array <double,ntuple> data164={40,8.4,10,0,52.5};dataID.push_back(data164);
      array <double,ntuple> data165={40,8.4,10,1,52.5};dataID.push_back(data165);
      array <double,ntuple> data166={40,8.4,15,0,45.5};dataID.push_back(data166);
      array <double,ntuple> data167={40,8.4,15,1,45.5};dataID.push_back(data167);
*/
      
      double n[dimNoise]={0.};
      double s[dimSignal]={0.};
      minimizeNoise(0,n,s);
      cout<<"For points: noise: ";
      for(int i=0;i<dimNoise;i++) cout<<n[i]<<", ";
      cout<<" signal: ";
      for(int i=0;i<dimSignal;i++) cout<<s[i]<<", ";
      cout<<endl;
      
      calculateRMS(n,s,0,true);
    }
    else{
      //double n[dimNoise]={10.2,0.99,-0.0086};
      //double s[dimSignal]={0,0,1};
      //double n[dimNoise]={13.6,1.0,-0.0013};
      //double n[dimNoise]={10.1,1.0,-0.0080};
      //double n[dimNoise]={9.9,1.0,-0.0080};
      //double s[dimSignal]={0,0,1};
      //double n[dimNoise]={13.77,1.00,-0.0034};
      //double s[dimSignal]={0,0,1};
      //double n[dimNoise]={2.49488, 1.1933, -0.0206882};
      //double s[dimSignal]={2.28612, 3.03863, 12.7583, 0.480369, 0.0057945};
      //double n[dimNoise]={141.475, -3.52019, 0.312033};
      double n[dimNoise]={970.062, 2.47523e-08, 1.77921};
      double s[dimSignal]={0.,0.,1.};
      //double n[dimNoise]={65.0199, 4.71792, 0.0260095, -0.00272133, 6.66942e-05};
      //double s[dimSignal]={0.,0.,1.};
      for(int ipc=0;ipc<nPC;ipc++){   
	for(int idr=0;idr<nDR;idr++){   
	  for(int ie=0;ie<nNRJ;ie++){
	    vertexResolution(event_energy[ie],vertex_resolution,efficiency,migrationProbability,nPMTtypes,monoNRJ,event_pc[ipc],event_dr[idr]);
	    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	      //cout<<"PMT type ="<<pmtType<<", Efficiency ="<<efficiency[pmtType]<<", vertex Resolution ="<<vertex_resolution[pmtType]<<endl;
	      vertex_resolution_energy[ipc][idr][pmtType][ie] = vertex_resolution[pmtType];
	      efficiency_energy[ipc][idr][pmtType][ie] = efficiency[pmtType];
	      migrationProbability_energy[ipc][idr][pmtType][ie] = migrationProbability[pmtType];
	      //signal_to_noise[ipc][idr][pmtType][ie]=SignalToNoise(event_energy[ie],event_pc[ipc],event_dr[idr],0,1,-0.01);
	      //signal_to_noise[ipc][idr][pmtType][ie]=SignalToNoise(event_energy[ie],event_pc[ipc],event_dr[idr],0,0.0076,-0.00021);
	      signal_to_noise[ipc][idr][pmtType][ie]=SignalToNoise(event_pc[ipc],event_dr[idr],event_energy[ie],n,s);
	      array <double,ntuple> data={event_pc[ipc],event_dr[ipc],event_energy[ie],pmtType,vertex_resolution[pmtType]};
	      dataID.push_back(data);
	      //cout<<"array <double,ntuple> data"<<counter<<"={"<<event_pc[ipc]<<","<<event_dr[ipc]<<","<<event_energy[ie]<<","<<pmtType<<","<<vertex_resolution[pmtType]<<"};dataID.push_back(data"<<counter<<");"<<endl;
	      counter++;
	    }
	  }
	}
      }
    }
    
    TGraph * genergy[nPC][nDR][nPMTtypes];
    TGraph * gefficiency[nPC][nDR][nPMTtypes];
    TGraph * gmigrationProbability[nPC][nDR][nPMTtypes];
    counter=0;
    //int color[12]={632, 800, 400, 416+2, 416, 416-8, 600, 616+2, 432, 1, 920, 616};
    //int color[nPC][nDR] = {{432,600,416,416+2},{800,632,616,880},{920,921,922,1}};
    //int color[nPC][nDR] = {{432,600},{800,632},{920,921}};
    //static const int ccolor[] = {632, 800, 400, 416+2, 416, 416-8, 600, 616+2, 432, 1, 920, 616};
    //std::vector <int> color (ccolor, ccolor + sizeof(ccolor) / sizeof(ccolor[0]));
    //int ic=0;
    //while(color.size() <NBinsTrueMom*NBinsTrueAngle){
    //color.push_back(921+ic);
    //ic++;
    //}
    
    TCanvas * cEnergy = new TCanvas("cEnergy","");
    TCanvas * cEnergy_perDR[nDR];
    TCanvas * cEnergy_perPC[nPC];
    for(int ipc=0;ipc<nPC;ipc++){
      for(int idr=0;idr<nDR;idr++){
	if(ipc==0) cEnergy_perDR[idr] = new TCanvas(Form("cEnergy_perDR_%d",idr),"");
	if(idr==0) cEnergy_perPC[ipc] = new TCanvas(Form("cEnergy_perPC_%d",ipc),"");
	//for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	  for(int pmtType=0;pmtType<1;pmtType++){
	    genergy[ipc][idr][pmtType] = new TGraph(nNRJ,event_energy,vertex_resolution_energy[ipc][idr][pmtType]);
	    genergy[ipc][idr][pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	    genergy[ipc][idr][pmtType]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	    genergy[ipc][idr][pmtType]->GetXaxis()->SetRangeUser(3,15);
	    genergy[ipc][idr][pmtType]->GetYaxis()->SetRangeUser(0,300);
	    genergy[ipc][idr][pmtType]->SetMarkerStyle(20);
	    genergy[ipc][idr][pmtType]->SetLineColor(color[ipc]);
	    genergy[ipc][idr][pmtType]->SetMarkerColor(color[ipc]);
	    genergy[ipc][idr][pmtType]->SetLineStyle(style[idr]);
	    cEnergy->cd();
	    if(counter==0) genergy[ipc][idr][pmtType]->Draw("ALP");
	    else genergy[ipc][idr][pmtType]->Draw("sameLP");
	    cEnergy_perDR[idr]->cd();
	    if(ipc==0) genergy[ipc][idr][pmtType]->Draw("ALP");
	    else genergy[ipc][idr][pmtType]->Draw("sameLP");
	    cEnergy_perPC[ipc]->cd();
	    if(idr==0) genergy[ipc][idr][pmtType]->Draw("ALP");
	    else genergy[ipc][idr][pmtType]->Draw("sameLP");
	  //genergy[ipc][idr][pmtType]->SaveAs(Form("plots/genergy[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	    counter++;
	  }
      }
    }
    cEnergy->SaveAs("plots/vertexEnergy.eps");
    for(int ipc=0;ipc<nPC;ipc++) cEnergy_perPC[ipc]->SaveAs(Form("plots/vertexEnergy_PC%d.eps",ipc));
    for(int idr=0;idr<nDR;idr++) cEnergy_perDR[idr]->SaveAs(Form("plots/vertexEnergy_DR%d.eps",idr));

    TCanvas * cSignalToNoise = new TCanvas("cSignalToNoise","");
    TGraph * gsignal_to_noise[nPC][nDR][nPMTtypes];
    counter=0;

    for(int ipc=0;ipc<nPC;ipc++){   
      for(int idr=0;idr<nDR;idr++){
	//for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	  for(int pmtType=0;pmtType<1;pmtType++){
	    gsignal_to_noise[ipc][idr][pmtType] = new TGraph(nNRJ,signal_to_noise[ipc][idr][pmtType],vertex_resolution_energy[ipc][idr][pmtType]);
	    gsignal_to_noise[ipc][idr][pmtType]->GetXaxis()->SetTitle("S/#sqrt(N)");
	    gsignal_to_noise[ipc][idr][pmtType]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	    //gsignal_to_noise[ipc][idr][pmtType]->GetXaxis()->SetRangeUser(3,15);
	    gsignal_to_noise[ipc][idr][pmtType]->GetYaxis()->SetRangeUser(0,300);
	    gsignal_to_noise[ipc][idr][pmtType]->SetMarkerStyle(20);
	    gsignal_to_noise[ipc][idr][pmtType]->SetLineColor(color[ipc]);
	    gsignal_to_noise[ipc][idr][pmtType]->SetMarkerColor(color[ipc]);
	    gsignal_to_noise[ipc][idr][pmtType]->SetLineStyle(style[idr]);
	    if(counter==0) gsignal_to_noise[ipc][idr][pmtType]->Draw("ALP");
	    else gsignal_to_noise[ipc][idr][pmtType]->Draw("sameLP");
	  //gsignal_to_noise[ipc][idr][pmtType]->SaveAs(Form("plots/gsignal_to_noise[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	  counter++;
	  }
      }
    }
    cSignalToNoise->SaveAs("plots/signalToNoise.eps");

    
    TCanvas * cEfficiency = new TCanvas("cEfficiency","");
    TCanvas * cEfficiency_perDR[nDR];
    TCanvas * cEfficiency_perPC[nPC];
    counter=0;
    for(int ipc=0;ipc<nPC;ipc++){
      for(int idr=0;idr<nDR;idr++){
	if(ipc==0) cEfficiency_perDR[idr] = new TCanvas(Form("cEfficiency_perDR_%d",idr),"");
	if(idr==0) cEfficiency_perPC[ipc] = new TCanvas(Form("cEfficiency_perPC_%d",ipc),"");
	//for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	  for(int pmtType=0;pmtType<1;pmtType++){
	    gefficiency[ipc][idr][pmtType] = new TGraph(nNRJ,event_energy,efficiency_energy[ipc][idr][pmtType]);
	    gefficiency[ipc][idr][pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	    gefficiency[ipc][idr][pmtType]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	    gefficiency[ipc][idr][pmtType]->GetXaxis()->SetRangeUser(3,15);
	    gefficiency[ipc][idr][pmtType]->GetYaxis()->SetRangeUser(0,1);
	    gefficiency[ipc][idr][pmtType]->SetMarkerStyle(20);
	    gefficiency[ipc][idr][pmtType]->SetLineColor(color[ipc]);
	    gefficiency[ipc][idr][pmtType]->SetMarkerColor(color[ipc]);
	    gefficiency[ipc][idr][pmtType]->SetLineStyle(style[idr]);
	    cEfficiency->cd();
	    if(counter==0) gefficiency[ipc][idr][pmtType]->Draw("ALP");
	    else gefficiency[ipc][idr][pmtType]->Draw("sameLP");
	    cEfficiency_perDR[idr]->cd();
	    if(ipc==0) gefficiency[ipc][idr][pmtType]->Draw("ALP");
	    else gefficiency[ipc][idr][pmtType]->Draw("sameLP");
	    cEfficiency_perPC[ipc]->cd();
	    if(idr==0) gefficiency[ipc][idr][pmtType]->Draw("ALP");
	    else gefficiency[ipc][idr][pmtType]->Draw("sameLP");
	  //gefficiency[ipc][idr][pmtType]->SaveAs(Form("plots/gefficiency[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	    counter++;
	  }
      }
    }
    cEfficiency->SaveAs("plots/vertexEfficiency.eps");
    for(int ipc=0;ipc<nPC;ipc++) cEfficiency_perPC[ipc]->SaveAs(Form("plots/vertexEfficiency_PC%d.eps",ipc));
    for(int idr=0;idr<nDR;idr++) cEfficiency_perDR[idr]->SaveAs(Form("plots/vertexEfficiency_DR%d.eps",idr));

      
    TCanvas * cMigrationProbability = new TCanvas("cMigrationProbability","");
    TCanvas * cMigrationProbability_perDR[nDR];
    TCanvas * cMigrationProbability_perPC[nPC];
    counter=0;
    for(int ipc=0;ipc<nPC;ipc++){
      for(int idr=0;idr<nDR;idr++){
	if(ipc==0) cMigrationProbability_perDR[idr] = new TCanvas(Form("cMigrationProbability_perDR_%d",idr),"");
	if(idr==0) cMigrationProbability_perPC[ipc] = new TCanvas(Form("cMigrationProbability_perPC_%d",ipc),"");
	//for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
	  for(int pmtType=0;pmtType<1;pmtType++){
	    gmigrationProbability[ipc][idr][pmtType] = new TGraph(nNRJ,event_energy,migrationProbability_energy[ipc][idr][pmtType]);
	    gmigrationProbability[ipc][idr][pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	    gmigrationProbability[ipc][idr][pmtType]->GetYaxis()->SetTitle("Migration probability");
	    gmigrationProbability[ipc][idr][pmtType]->GetXaxis()->SetRangeUser(3,15);
	    gmigrationProbability[ipc][idr][pmtType]->GetYaxis()->SetRangeUser(0,0.4);
	    gmigrationProbability[ipc][idr][pmtType]->SetMarkerStyle(20);
	    gmigrationProbability[ipc][idr][pmtType]->SetLineColor(color[ipc]);
	    gmigrationProbability[ipc][idr][pmtType]->SetMarkerColor(color[ipc]);
	    gmigrationProbability[ipc][idr][pmtType]->SetLineStyle(style[idr]);
	    cMigrationProbability->cd();
	    if(counter==0) gmigrationProbability[ipc][idr][pmtType]->Draw("ALP");
	    else gmigrationProbability[ipc][idr][pmtType]->Draw("sameLP");
	    cMigrationProbability_perDR[idr]->cd();
	    if(ipc==0) gmigrationProbability[ipc][idr][pmtType]->Draw("ALP");
	    else gmigrationProbability[ipc][idr][pmtType]->Draw("sameLP");
	    cMigrationProbability_perPC[ipc]->cd();
	    if(idr==0) gmigrationProbability[ipc][idr][pmtType]->Draw("ALP");
	    else gmigrationProbability[ipc][idr][pmtType]->Draw("sameLP");
	  //gmigrationProbability[ipc][idr][pmtType]->SaveAs(Form("plots/gmigrationProbability[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	    counter++;
	  }
      }
    }
    cMigrationProbability->SaveAs("plots/vertexMigrationProbability.eps");
    for(int ipc=0;ipc<nPC;ipc++) cMigrationProbability_perPC[ipc]->SaveAs(Form("plots/vertexMigrationProbability_PC%d.eps",ipc));
    for(int idr=0;idr<nDR;idr++) cMigrationProbability_perDR[idr]->SaveAs(Form("plots/vertexMigrationProbability_DR%d.eps",idr));

    /*    
    TCanvas * cEfficiency = new TCanvas("cEfficiency","");
    TGraph * gefficiency[nPMTtypes];
    for(int pmtType=0;pmtType<nPMTtypes;pmtType++){
      gefficiency[pmtType] = new TGraph(nNRJ,event_energy,efficiency_energy[pmtType]);
      gefficiency[pmtType]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
      gefficiency[pmtType]->GetYaxis()->SetTitle("Signal efficiency");
      gefficiency[pmtType]->GetXaxis()->SetRangeUser(3,15);
      gefficiency[pmtType]->GetYaxis()->SetRangeUser(0,1.05);
      gefficiency[pmtType]->SetMarkerStyle(20);
      if(pmtType==0){
	gefficiency[pmtType]->SetLineColor(kBlue);
	gefficiency[pmtType]->SetMarkerColor(kBlue);
	gefficiency[pmtType]->Draw("ALP");
      }
      else{
	gefficiency[pmtType]->SetLineColor(kRed);
	gefficiency[pmtType]->SetMarkerColor(kRed);
	gefficiency[pmtType]->Draw("sameLP");
      }
      gefficiency[pmtType]->SaveAs(Form("plots/gefficiency_pmtType%d_%s.root",pmtType,filename));
    }
    cEfficiency->SaveAs("plots/vertexEfficiency.eps");    

    */
  }
  
  theApp.Run();
  return 0;
}

  
