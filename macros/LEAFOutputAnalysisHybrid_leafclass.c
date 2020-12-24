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


//Assumes that function is of the form (a + b*sqrt(N) + c*sqrt(N)^2)/S
double calculateRMS(double n[dimNoise], double s[dimSignal], bool drawOption=false){

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
    double vertex_resolution=data[3];
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
  //double x[npoints]={0.};
  //double ex[npoints]={0.};
  //double y[npoints]={0.};
  //double ey[npoints]={0.};
  double max_signal_to_noise=-1e9;
  double min_signal_to_noise=1e9;

  int count=0;
  for(int i=0;i<dataID.size();i++){
    array < double,ntuple > data=dataID[i];
    double pc=data[0];
    double dr=data[1];
    double e=data[2];
    double vertex_resolution=data[3];
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

  RMS=calculateRMS(noise,signal);//

}



//here is the function where minimization gonna take place. The loop is inside this function
void minimizeNoise(double n[dimNoise], double s[dimSignal]){
  
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

  for(int i=0;i<dimSignal;i++) minimizer->FixParameter(dimNoise+i);
  Mig[1]=1e-8;minimizer->ExecuteCommand("SIMPLEX",Mig,2);
  Mig[1]=1e-1;minimizer->ExecuteCommand("MIGRAD",Mig,2);

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

int vertexResolution(double energy, double vertex_resolution, double efficiency, double misIDFV, bool monoEnergy=true, double photocov=40., double darkrate=4.2){


  
  //##############################################################################################
  //##############################################################################################
  //############################# BASIC CONSTANTS AND CUT VALUES #################################
  //##############################################################################################
  //##############################################################################################
  //Tank constants
  double tankRadius = 7080/2;//In cm
  double tankHeight = 5480/2;//In cm
  //Cuts
  bool FVCut = true;
  double FVCutValue = 0;//in cm
  bool ExtCut = false;
  double ExtCutvalue = 1200;//in cm  
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################

  

  
  //##############################################################################################
  //##############################################################################################
  //##################################### OPENING THE INPUT FILE #################################
  //##############################################################################################
  //##############################################################################################
  bool verbose=true;
  //TFile *file = new TFile(Form("LEAFv2_hk20pc_4.2kHzBL_fulltank_nodirectionality_%dMeV_v2.root",((int) energy) ));
  TFile *file = new TFile(Form("LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzBL_fulltank_nodirectionality_%dMeV_v2.root",((int) energy) ));
  //TFile *file = new TFile(Form("LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzBL_closewall_nodirectionality_%dMeV_v2.root",((int) energy) ));
  if(verbose) cout<<"We will try to open file = "<<file->GetName()<<endl;
  file->cd();
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################



  
  //##############################################################################################
  //##############################################################################################
  //#################### INPUT TTREE LINKING TO VARIABLES ########################################
  //##############################################################################################
  //##############################################################################################
  TTree *bstree = (TTree*)file->Get("Reduced");
  int nentries = bstree->GetEntries();
  if(verbose) cout<<"Number of entries = "<<nentries<<endl;

  //Preparation of host variables
  std::vector<float> *true_vertex_X = 0;std::vector<float> *true_vertex_Y = 0;std::vector<float> *true_vertex_Z = 0;std::vector<float> *true_vertex_T = 0;
  std::vector<float> *true_dir_X = 0;std::vector<float> *true_dir_Y = 0;std::vector<float> *true_dir_Z = 0;
  double lf_time = 0.;
  double lf_vertex[3] = {0.};
  double lf_dir[3][3] = {{0.}};
  double lf_wall = 0;
  int lf_n50[3] = {0};
  double lf_good = 0.;
  double lf_dirKS[3] = {0.};
  int lf_intime = 0;
  int eventId = 0, triggerId = 0, ID_hits_200 = 0, ID_hits = 0;
  
  //True variables
  bstree->SetBranchAddress("true_vertex_X", &true_vertex_X);
  bstree->SetBranchAddress("true_vertex_Y", &true_vertex_Y);
  bstree->SetBranchAddress("true_vertex_Z", &true_vertex_Z);
  bstree->SetBranchAddress("true_vertex_T", &true_vertex_T);
  bstree->SetBranchAddress("true_dir_X", &true_dir_X);
  bstree->SetBranchAddress("true_dir_Y", &true_dir_Y);
  bstree->SetBranchAddress("true_dir_Z", &true_dir_Z);  
  //LEAF reconstructed variables			     
  bstree->SetBranchAddress("lf_time", &lf_time);
  bstree->SetBranchAddress("lf_vertex",lf_vertex);
  bstree->SetBranchAddress("lf_dir",lf_dir);
  bstree->SetBranchAddress("lf_good", &lf_good);
  bstree->SetBranchAddress("lf_dirKS", lf_dirKS);
  bstree->SetBranchAddress("lf_wall",&lf_wall);
  bstree->SetBranchAddress("lf_n50", lf_n50);
  //Pure detector raw variables
  bstree->SetBranchAddress("eventId",&eventId);
  bstree->SetBranchAddress("triggerId",&triggerId);
  bstree->SetBranchAddress("ID_hits", &ID_hits);
  bstree->SetBranchAddress("ID_hits_200",&ID_hits_200);
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################





  
  //##############################################################################################
  //##############################################################################################
  //############################# DECLARATION OF OUTPUT HISTOGRAMS  ##############################
  //##############################################################################################
  //##############################################################################################

  //############################## Histo declaration #############################################
  //Vertex related
  TH2D * hvertexLongitudinalXTime;
  TH1D * hvertexLongitudinal;TH1D * hvertexOrtho;
  TH1D * hvertex;
  TH1D * hradius;
  double nVertexFound = 0., nEvents = 0., nEventsTotal = 0.;
  TH2D * hvertex2D;
  TH2D * htruevertex2D;
  TH2D * hvertexXdwall;
  TH2D * hvertexXtowall;
  TH2D * hdwallXdwallrec; 

  //Nhits related
  TH1D * hnhits;
  TH2D * hdwallXnhits;
  TH2D * htowallXnhits;

  //Direction related
  TH1D * hdirection;
  
  //Goodness of fit & dirKS related
  TH1D * hgood;
  TH1D * hdirKS;
  
  //Mix of previous variables
  TH2D * hvertexXnhits_intime;
  TH2D * hvertexXgood;
  TH2D * hvertexXnhits;
 
  //Cut related
  TH1D * hvertexMisIDFV; TH1D * hvertexMisIDFV_total;
  TH1D * hvertexMisIDExtCut; TH1D * hvertexMisIDExtCut_total;
  double nMisIDExtCut = 0., nMisIDExtCut_total = 0.;
  double nMisIDFV = 0., nMisIDFV_total = 0.;
  //###########################################################################################

  

  //############################## Histo binning choice #######################################
  //Binning chosen
  int nbin_vertex = 6000;//Number of bins used for vertex resolution
  double max_vertex = 6000;
  int nbin_direction = 180;//Number of bins used for angle resolution
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
  int nbin_good = 100;//Number of bins used for vertex resolution
  double max_good = 1000;
  int nbin_dirKS = 50;//Number of bins used for vertex resolution
  double max_dirKS = 1;
  //###########################################################################################


  //############################## Histo initialization #######################################

  //Vertex related
  hvertexLongitudinalXTime = new TH2D(Form("hvertexLongitudinalXTime"),"",nbin_vertex,0,max_vertex,nbin_time,-max_time,max_time);
  hvertexLongitudinalXTime->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexLongitudinalXTime->GetYaxis()->SetTitle("dwall (cm)");			     

  hvertexLongitudinal = new TH1D(Form("hvertexLongitudinal"),"",nbin_vertex,0,max_vertex);
  hvertexLongitudinal->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hvertexOrtho = new TH1D(Form("hvertexOrtho"),"",nbin_vertex,0,max_vertex);
  hvertexOrtho->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hvertex = new TH1D(Form("hvertex"),"",nbin_vertex,0,max_vertex);
  hvertex->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");

  hradius = new TH1D(Form("hradius"),"",nbin_vertex,0,max_vertex*max_vertex);
  hradius->GetXaxis()->SetTitle("R2 (cm2)");
 
  hvertex2D = new TH2D(Form("hvertex2D"),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
  hvertex2D->GetXaxis()->SetTitle("x position (cm)");
  hvertex2D->GetYaxis()->SetTitle("y position (cm)");

  htruevertex2D = new TH2D(Form("htruevertex2D"),"",nbin_tanksize,min_tanksize,max_tanksize,nbin_tanksize,min_tanksize,max_tanksize);
  htruevertex2D->GetXaxis()->SetTitle("x position (cm)");
  htruevertex2D->GetYaxis()->SetTitle("y position (cm)");

  hdwallXdwallrec = new TH2D(Form("hdwallXdwallrec"),"",nbin_distance,0,max_distance,20,-1,1);
  hdwallXdwallrec->GetXaxis()->SetTitle("dwall true (cm)");
  hdwallXdwallrec->GetYaxis()->SetTitle("dwall rec  - dwall true / dwall true");

  hvertexXdwall = new TH2D(Form("hvertexXdwall"),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
  hvertexXdwall->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXdwall->GetYaxis()->SetTitle("dwall (cm)");

  hvertexXtowall = new TH2D(Form("hvertexXtowall"),"",nbin_vertex,0,max_vertex,nbin_distance,0,max_distance);
  hvertexXtowall->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXtowall->GetYaxis()->SetTitle("towall (cm)");

  
  //Nhits related
  hnhits = new TH1D(Form("hnhits"),"",nbin_nhits,0,max_nhits);
  hnhits->GetXaxis()->SetTitle("number of hits");

  hdwallXnhits = new TH2D(Form("hdwallXnhits"),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
  hdwallXnhits->GetXaxis()->SetTitle("dwall (cm)");
  hdwallXnhits->GetYaxis()->SetTitle("number of hits");

  htowallXnhits = new TH2D(Form("htowallXnhits"),"",nbin_distance,0,max_distance,nbin_nhits,0,max_nhits);
  htowallXnhits->GetXaxis()->SetTitle("towall (cm)");
  htowallXnhits->GetYaxis()->SetTitle("number of hits");

  
  //Direction related
  hdirection = new TH1D(Form("hdirection"),"",nbin_direction,0,max_direction);
  hdirection->GetXaxis()->SetTitle("angle from rec. to mc. direction (#circ)");


  //Goodness of fit & dirKS related
  hgood = new TH1D(Form("hgood"),"",nbin_good,0,max_good);
  hgood->GetXaxis()->SetTitle("good of fit");
  hdirKS = new TH1D(Form("hdirKS"),"",nbin_dirKS,0,max_dirKS);
  hdirKS->GetXaxis()->SetTitle("dirKS");

  
  //Mix of previous variables
  hvertexXnhits_intime = new TH2D(Form("hvertexXnhits_intime"),"",nbin_vertex,0,max_vertex,nbin_nhits,0.,max_nhits);
  hvertexXnhits_intime->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXnhits_intime->GetYaxis()->SetTitle("number of hits");

  hvertexXgood = new TH2D(Form("hvertexXgoddness"),"",nbin_vertex,0,max_vertex,nbin_good,0.,max_good);
  hvertexXgood->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXgood->GetYaxis()->SetTitle("good of fit");

  hvertexXnhits = new TH2D(Form("hvertexXnhits"),"",nbin_vertex,0,max_vertex,nbin_nhits,0,max_nhits);
  hvertexXnhits->GetXaxis()->SetTitle("distance from rec. to mc. vertex (cm)");
  hvertexXnhits->GetYaxis()->SetTitle("number of hits");

  
  //Cut related
  hvertexMisIDFV = new TH1D(Form("hvertexMisIDFV"),"",20,0,200);
  hvertexMisIDFV->GetXaxis()->SetTitle("FV distance to the wall");
  hvertexMisIDFV->GetYaxis()->SetTitle("Probability of migration from out to in FV");
  hvertexMisIDFV_total = new TH1D(Form("hvertexMisIDFV_total"),"",20,0,200);

  hvertexMisIDExtCut = new TH1D(Form("hvertexMisIDExtCut"),"",20,0,2000);
  hvertexMisIDExtCut->GetXaxis()->SetTitle("ExtCut distance to the wall");
  hvertexMisIDExtCut->GetYaxis()->SetTitle("Probability of migration from out to in ExtCut");
  hvertexMisIDExtCut_total = new TH1D(Form("hvertexMisIDExtCut_total"),"",20,0,2000);
 
  //###########################################################################################

  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################
  //##############################################################################################



  
 




  
  //######################################################################################################
  //######################################################################################################
  //############################## LOOP OVER EACH DAQ/TRUE VERTEX EVENT ################################## 
  //######################################################################################################
  //######################################################################################################


  for (int entry=0; entry<nentries; entry++) {
    bstree->GetEntry(entry);
    nEventsTotal++;
    
    //######################################################################################################
    //######################################################################################################
    //###################################### LOOP OVER EACH TRIGGERS #######################################
    //######################################################################################################
    //######################################################################################################
    
    //######################################################################################################
    // One DAQ/Sim event can contain several vertices, and several trigger per vertices. Here, we assume one
    // vertex per event, but allow several trigger for each of them.
    //######################################################################################################
    int ntriggers = 1;//true_vertex_X->size();
    if(ntriggers!=1) cout<<"In event "<<eventId<<", there are "<<ntriggers<<" triggers"<<endl;
    for(int itr=0;itr<ntriggers;itr++){
      float mc_vertex[3] = {true_vertex_X->at(itr),true_vertex_Y->at(itr),true_vertex_Z->at(itr)};
      float mc_time = true_vertex_T->at(itr);
      float mc_dir[3] = {true_dir_X->at(itr),true_dir_Y->at(itr),true_dir_Z->at(itr)};
      if(verbose && entry%10000==0){
	std::cout<<"####################################"<<std::endl<<"############ True event ############"<<std::endl;
	std::cout<<"Vertex pos = "<<mc_vertex[0]<<", "<<mc_vertex[1]<<", "<<mc_vertex[2]<<endl;
	std::cout<<"Particle dir = "<<mc_dir[0]<<", "<<mc_dir[1]<<", "<<mc_dir[2]<<endl;
	std::cout<<"############ Reco event ############"<<std::endl;
	std::cout<<"Vertex pos = "<<lf_vertex[0]<<", "<<lf_vertex[1]<<", "<<lf_vertex[2]<<endl;
	//std::cout<<"Particle dir = "<<lf_dir[2][0]<<", "<<lf_dir[2][1]<<", "<<lf_dir[2][2]<<endl;
	std::cout<<"Particle dir PMT 1 = "<<lf_dir[0][0]<<", "<<lf_dir[0][1]<<", "<<lf_dir[0][2]<<endl;
	std::cout<<"Particle dir PMT 2 = "<<lf_dir[1][0]<<", "<<lf_dir[1][1]<<", "<<lf_dir[1][2]<<endl;
	std::cout<<"Particle dir PMT 3 = "<<lf_dir[2][0]<<", "<<lf_dir[2][1]<<", "<<lf_dir[2][2]<<endl;
	std::cout<<"####################################"<<std::endl<<std::endl;
      }
      //######################################################################################################
      //######################################################################################################
      //############################# Calculation of dWall and to wall #######################################
      //######################################################################################################
      //######################################################################################################
      
      double radius_true_vertex = TMath::Sqrt( pow(mc_vertex[0],2) + pow(mc_vertex[1],2) );
      //cout<<"Radius = "<<TMath::Sqrt( pow(mc_vertex[0],2) + pow(mc_vertex[1],2) )<<", and z position = "<<mc_vertex[2]<<endl;
      //dwall is the smallest of distance from vertex to top/bottom vs to edge
      double dwall = TMath::Min(tankHeight-TMath::Abs(mc_vertex[2]),TMath::Abs(radius_true_vertex-tankRadius));
      htruevertex2D->Fill(mc_vertex[0],mc_vertex[1]);
      hvertex2D->Fill(lf_vertex[0],lf_vertex[1]);
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

      if(FVCut && lf_wall<FVCutValue) continue;
      //if(FVCut && dwall<FVCutValue) continue;
      if(ExtCut && pmt_to_vertex<ExtCutvalue) continue;
      nEvents++;

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
      hvertex->Fill(vertex_distance);
      hradius->Fill(pow(radius_rec_vertex,2));
      hnhits->Fill(ID_hits);
      hvertexXnhits->Fill(vertex_distance,ID_hits);
      hdwallXnhits->Fill(dwall,ID_hits);
      htowallXnhits->Fill(towall,ID_hits);
      hvertexXnhits_intime->Fill(vertex_distance,lf_intime);
      nVertexFound++;

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
      hvertexLongitudinal->Fill(TMath::Abs(vector_director_longitudinal));
      hvertexOrtho->Fill(vector_director_ortho);
      hvertexLongitudinalXTime->Fill(TMath::Abs(vector_director_longitudinal),lf_vertex[3]-mc_vertex[3]);

      //######################################################################################################
      // Vertex resolution vs dwall and towall
      //######################################################################################################
      hvertexXdwall->Fill(vertex_distance,dwall);
      hvertexXtowall->Fill(vertex_distance,towall);
      hdwallXdwallrec->Fill(dwall,(dwall_rec-dwall)/dwall);
     
      for(int ibinx=1;ibinx<=hvertexMisIDFV->GetNbinsX();ibinx++){
	double dwall_FV_definition=hvertexMisIDFV->GetBinLowEdge(ibinx)+hvertexMisIDFV->GetBinWidth(ibinx);
	if(dwall<dwall_FV_definition){//check only if the event is out of FV
	  bool misID=migrateFromOVtoFV(dwall, dwall_rec, FVCutValue);
	  if(misID) hvertexMisIDFV->Fill(dwall_FV_definition-hvertexMisIDFV->GetBinWidth(ibinx));
	  hvertexMisIDFV_total->Fill(dwall_FV_definition-hvertexMisIDFV->GetBinWidth(ibinx));
	}
      }
    
      if(dwall<FVCutValue){
	bool migration=migrateFromOVtoFV(dwall,dwall_rec,FVCutValue);
	if(migration) nMisIDFV++;
	nMisIDFV_total++;
      }

      for(int ibinx=1;ibinx<=hvertexMisIDExtCut->GetNbinsX();ibinx++){
	double ExtCut_definition=hvertexMisIDExtCut->GetBinLowEdge(ibinx)+hvertexMisIDExtCut->GetBinWidth(ibinx);
	//cout<<dwall_ExtCut_definition<<endl;
	//if(pmt_to_vertex<ExtCut_definition){//check only if the event is out of ExtCut
	bool misID=migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_rec, ExtCut_definition);
	if(misID) hvertexMisIDExtCut->Fill(ExtCut_definition-hvertexMisIDExtCut->GetBinWidth(ibinx));
	hvertexMisIDExtCut_total->Fill(ExtCut_definition-hvertexMisIDExtCut->GetBinWidth(ibinx));
	//cout<<"misID="<<misID<<endl;
	//}
      }

      //if(pmt_to_vertex<ExtCutvalue){
      bool migration=migrateFromOVtoFV(pmt_to_vertex, pmt_to_vertex_rec, ExtCutvalue);
      //bool migration=(pmt_to_vertex>ExtCutvalue);
      if(migration) nMisIDExtCut++;
      nMisIDExtCut_total++;

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
      double direction_angle = TMath::ACos(lf_dir[2][0]*mc_dir[0]+lf_dir[2][1]*mc_dir[1]+lf_dir[2][2]*mc_dir[2])*180/TMath::Pi();//Scalar product between true and rec dir.
      //if(direction_angle!=90.)
      hdirection->Fill(direction_angle);
      
      if(verbose && direction_angle==90.){
	std::cout<<"####################################"<<std::endl<<"Direction difference rec vs true = "<<direction_angle<<"degrees"<<std::endl<<"############ True event ############"<<std::endl;
	std::cout<<"Vertex pos = "<<mc_vertex[0]<<", "<<mc_vertex[1]<<", "<<mc_vertex[2]<<endl;
	std::cout<<"Particle dir = "<<mc_dir[0]<<", "<<mc_dir[1]<<", "<<mc_dir[2]<<endl;
	std::cout<<"############ Reco event ############"<<std::endl;
	std::cout<<"Vertex pos = "<<lf_vertex[0]<<", "<<lf_vertex[1]<<", "<<lf_vertex[2]<<endl;
	std::cout<<"Particle dir PMT 1 = "<<lf_dir[0][0]<<", "<<lf_dir[0][1]<<", "<<lf_dir[0][2]<<endl;
	std::cout<<"Particle dir PMT 2 = "<<lf_dir[1][0]<<", "<<lf_dir[1][1]<<", "<<lf_dir[1][2]<<endl;
	std::cout<<"Particle dir PMT 3 = "<<lf_dir[2][0]<<", "<<lf_dir[2][1]<<", "<<lf_dir[2][2]<<endl;
	std::cout<<"dwall = "<<lf_wall<<" va rec dwall = "<<dwall_rec<<", n50 hits = "<<lf_n50[2]<<std::endl;
	std::cout<<"####################################"<<std::endl<<std::endl;
      }
      
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################
      //######################################################################################################

      //######################################################################################################
      //######################################################################################################
      //################################# Goodness of fit & dirKS ############################################
      //######################################################################################################
      //######################################################################################################

      //######################################################################################################
      // Simple goodness of fit
      //######################################################################################################
      hgood->Fill(lf_good);

      //######################################################################################################
      // Simple dirKS of fit
      //######################################################################################################
      hdirKS->Fill(lf_dirKS[2]);
      
      
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



  
  efficiency = (nEvents==0?0:nVertexFound/nEvents);
  misIDFV = (nEvents==0?0:nMisIDFV/nMisIDFV_total);
  cout << "Efficiency in loop : "<< efficiency << " and misID = " << misIDFV << endl;

  double misIDExtCut;
  misIDExtCut = (nEvents==0?0:nMisIDExtCut/nMisIDExtCut_total);
  cout<<"Efficiency in loop : "<< efficiency << " and misID = "<<misIDExtCut<<endl;
  cout<<"Number of events passing the FV and other cuts = "<<nEvents/nEventsTotal<<endl;

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  vertex_resolution = vRes(hvertex);
  double error = vResAndError(hvertex,CredibleInterval);
  cout<<"Vertex resolution = "<<vertex_resolution<<endl;

  TCanvas * cvertex2D = new TCanvas("cvertex2D","");
  hvertex2D->Draw("colz");
  cvertex2D->SaveAs("plots/vertex2D.eps");
  

  TCanvas * ctruevertex2D = new TCanvas("ctruevertex2D","");
  htruevertex2D->Draw("colz");
  ctruevertex2D->SaveAs("plots/truevertex2D.eps");

  
  TCanvas * cvertex = new TCanvas("cvertex","");
  hvertex->Rebin(10);
  hvertex->GetXaxis()->SetRangeUser(0,300);
  hvertex->Scale(1./hvertex->Integral());
  hvertex->SetLineColor(kBlue);
  hvertex->Draw();
  hvertex->SaveAs(Form("root_output/hvertex_%s.root",file->GetName()));
  cvertex->SaveAs("plots/vertex.eps");

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexLongitudinal = new TCanvas("cvertexLongitudinal","");
  hvertexLongitudinal->Rebin(10);
  hvertexLongitudinal->GetXaxis()->SetRangeUser(0,300);
  hvertexLongitudinal->SetLineColor(kBlue);
  hvertexLongitudinal->Draw();
  cvertexLongitudinal->SaveAs("plots/vertexLongitudinal.eps");
  
  TCanvas * cvertexOrtho = new TCanvas("cvertexOrtho","");
  hvertexOrtho->Rebin(10);
  hvertexOrtho->GetXaxis()->SetRangeUser(0,300);
  hvertexOrtho->SetLineColor(kBlue);
  hvertexOrtho->Draw();
  cvertexOrtho->SaveAs("plots/vertexOrtho.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexLongitudinalXTime = new TCanvas("cvertexLongitudinalXTime","");
  hvertexLongitudinalXTime->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////



  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall = new TCanvas("cvertexXdwall","");
  hvertexXdwall->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall = new TCanvas("cvertexXtowall","");
  hvertexXtowall->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall1DAll = new TCanvas("cvertexXdwall1DAll","");
  const int nSlicesdwall = nbin_distance;//16;//hvertexXdwall[0]->GetNbinsX();
  cvertexXdwall1DAll->Divide(TMath::Sqrt(nSlicesdwall),TMath::Sqrt(nSlicesdwall));
  TH1D * hvertexXdwall1D[nSlicesdwall];
  double dwallPosition[nSlicesdwall];
  double resolution68[nSlicesdwall];
  double error_resolution68[nSlicesdwall];
  double error_x[nSlicesdwall];
  TGraphErrors * gdwall;
  for(int ibinx=1;ibinx<=nSlicesdwall;ibinx++){
    hvertexXdwall1D[ibinx-1] = (TH1D*) hvertexXdwall->ProjectionX(Form("hvertexXdwall1D_%d",ibinx-1),ibinx,ibinx);
    dwallPosition[ibinx-1]=hvertexXdwall->GetYaxis()->GetBinCenter(ibinx);
    resolution68[ibinx-1]=vRes(hvertexXdwall1D[ibinx-1]);
    error_resolution68[ibinx-1]=vResAndError(hvertexXdwall1D[ibinx-1],CredibleInterval);
    //cout<<"Error = "<<error_resolution68[ibinx-1]<<endl;
    error_x[ibinx-1]=0;
    //hvertexXdwall1D[ibinx-1]->Scale(1./hvertexXdwall1D[ibinx-1]->Integral());
    hvertexXdwall1D[ibinx-1]->GetXaxis()->SetTitle(hvertexXdwall->GetYaxis()->GetTitle());
    cvertexXdwall1DAll->cd(ibinx);
    hvertexXdwall1D[ibinx-1]->SetLineColor(kBlue);
    hvertexXdwall1D[ibinx-1]->Draw();
  }
  gdwall = new TGraphErrors(nSlicesdwall,dwallPosition,resolution68,error_x,error_resolution68);
  gdwall->GetXaxis()->SetTitle(hvertexXdwall->GetYaxis()->GetTitle());
  gdwall->GetYaxis()->SetTitle(hvertexXdwall->GetXaxis()->GetTitle());
  gdwall->SaveAs(Form("plots/gdwall_%s.root",file->GetName()));
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXdwall1D = new TCanvas("cvertexXdwall1D","");
  gdwall->SetMarkerStyle(20);
  gdwall->GetXaxis()->SetRangeUser(min_distance,max_distance);
  gdwall->GetYaxis()->SetRangeUser(50,100);
  gdwall->SetLineColor(kBlue);
  gdwall->SetMarkerColor(kBlue);
  gdwall->Draw("ALP");
  gdwall->SaveAs(Form("plots/gdwall_%s.root",file->GetName()));
  cvertexXdwall1D->SaveAs("plots/vertexXdwall.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////




  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall1DAll = new TCanvas("cvertexXtowall1DAll","");
  const int nSlicestowall = nbin_distance;//16;//hvertexXtowall[0]->GetNbinsX();
  cvertexXtowall1DAll->Divide(TMath::Sqrt(nSlicestowall),TMath::Sqrt(nSlicestowall));
  TH1D * hvertexXtowall1D[nSlicestowall];
  double towallPosition[nSlicestowall];
  double resolution68_towall[nSlicestowall];
  double error_resolution68_towall[nSlicestowall];
  double error_x_towall[nSlicestowall];
  TGraphErrors * gtowall;
  for(int ibinx=1;ibinx<=nSlicestowall;ibinx++){
    hvertexXtowall1D[ibinx-1] = (TH1D*) hvertexXtowall->ProjectionX(Form("hvertexXtowall1D_%d",ibinx-1),ibinx,ibinx);
    towallPosition[ibinx-1]=hvertexXtowall->GetYaxis()->GetBinCenter(ibinx);
    resolution68_towall[ibinx-1]=vRes(hvertexXtowall1D[ibinx-1]);
    error_resolution68_towall[ibinx-1]=vResAndError(hvertexXtowall1D[ibinx-1],CredibleInterval);
    //cout<<"Error = "<<error_resolution68_towall[ibinx-1]<<endl;
    error_x_towall[ibinx-1]=0;
    //hvertexXtowall1D[ibinx-1]->Scale(1./hvertexXtowall1D[ibinx-1]->Integral());
    hvertexXtowall1D[ibinx-1]->GetXaxis()->SetTitle(hvertexXtowall->GetYaxis()->GetTitle());
    cvertexXtowall1DAll->cd(ibinx);
    hvertexXtowall1D[ibinx-1]->SetLineColor(kBlue);
    hvertexXtowall1D[ibinx-1]->Draw();
    gtowall = new TGraphErrors(nSlicestowall,towallPosition,resolution68_towall,error_x_towall,error_resolution68_towall);
    gtowall->GetXaxis()->SetTitle(hvertexXtowall->GetYaxis()->GetTitle());
    gtowall->GetYaxis()->SetTitle(hvertexXtowall->GetXaxis()->GetTitle());
    gtowall->SaveAs(Form("plots/gtowall_%s.root",file->GetName()));
  }
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXtowall1D = new TCanvas("cvertexXtowall1D","");
  gtowall->SetMarkerStyle(20);
  gtowall->GetXaxis()->SetRangeUser(min_distance,max_distance);
  gtowall->GetYaxis()->SetRangeUser(50,100);
  gtowall->SetLineColor(kBlue);
  gtowall->SetMarkerColor(kBlue);
  gtowall->Draw("ALP");
  gtowall->SaveAs(Form("plots/gtowall_%s.root",file->GetName()));
  cvertexXtowall1D->SaveAs("plots/vertexXtowall.eps");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  TCanvas * cradius = new TCanvas("cradius","");
  hradius->Rebin(10);
  //hradius->GetXaxis()->SetRangeUser(0,300);
  hradius->Scale(1./hradius->Integral());
  hradius->SetLineColor(kBlue);
  hradius->Draw();
  hradius->SaveAs(Form("root_output/hradius_%s.root",file->GetName()));
  cradius->SaveAs("plots/radius.eps");


  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdwallXdwallrec = new TCanvas("cdwallXdwallrec","");
  //hdwallXdwallrec->Rebin2D(1,;
  hdwallXdwallrec->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cnhits = new TCanvas("cnhits","");
  hnhits->SetLineColor(kBlue);
  hnhits->Draw();
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXnhits = new TCanvas("cvertexXnhits","");
  hvertexXnhits->Rebin2D(1,2);
  hvertexXnhits->Draw("colz");  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdwallXnhits = new TCanvas("cdwallXnhits","");
  hdwallXnhits->Rebin2D(1,2);
  hdwallXnhits->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////

  
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * ctowallXnhits = new TCanvas("ctowallXnhits","");
  htowallXnhits->Rebin2D(1,2);
  htowallXnhits->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXnhits_intime = new TCanvas("cvertexXnhits_intime","");
  hvertexXnhits_intime->Rebin2D(10,1);
  hvertexXnhits_intime->GetYaxis()->SetRangeUser(0,10);
  hvertexXnhits_intime->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexXgood = new TCanvas("cvertexXgood","");
  hvertexXgood->Rebin2D(100,2);
  hvertexXgood->GetYaxis()->SetRangeUser(0.,1.);
  hvertexXgood->Draw("colz");
  ////////////////////////////////////////////////
  ////////////////////////////////////////////////


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdirection = new TCanvas("cdirection","");
  cout<<"Direction resolution = "<<vRes(hdirection)<<endl;
  //hdirection->Rebin(2);
  hdirection->GetXaxis()->SetRangeUser(0,100);
  hdirection->Scale(1./hdirection->Integral());
  hdirection->SetLineColor(kBlue);
  hdirection->Draw();
  hdirection->SaveAs(Form("plots/hdirection_%s.root",file->GetName()));
  cdirection->SaveAs("plots/direction.eps");

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cgood = new TCanvas("cgood","");
  //hgood->Rebin(2);
  hgood->Scale(1./hgood->Integral());
  hgood->SetLineColor(kBlue);
  hgood->Draw();
  hgood->SaveAs(Form("plots/hgood_%s.root",file->GetName()));
  cgood->SaveAs("plots/good.eps");

  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cdirKS = new TCanvas("cdirKS","");
  //hdirKS->Rebin(2);
  hdirKS->Scale(1./hdirKS->Integral());
  hdirKS->SetLineColor(kBlue);
  hdirKS->Draw();
  hdirKS->SaveAs(Form("plots/hdirKS_%s.root",file->GetName()));
  cdirKS->SaveAs("plots/dirKS.eps");


  ////////////////////////////////////////////////
  ////////////////////////////////////////////////
  TCanvas * cvertexMisIDFV = new TCanvas("cvertexMisIDFV","");
  gStyle->SetOptStat(0.);
  cvertexMisIDFV->SetGridx(true);
  cvertexMisIDFV->SetGridy(true);
  hvertexMisIDFV->Divide(hvertexMisIDFV_total);
  hvertexMisIDFV->SetLineWidth(2);
  hvertexMisIDFV->SetLineColor(kBlue);
  hvertexMisIDFV->GetXaxis()->SetRangeUser(0,200); 
  hvertexMisIDFV->GetYaxis()->SetRangeUser(0,1.); 
  hvertexMisIDFV->Draw();
  hvertexMisIDFV->SaveAs(Form("root_output/hvertexMisIDFV_%s.root",file->GetName()));
  cvertexMisIDFV->SaveAs("plots/vertexMisIDFV.eps");


  TCanvas * cvertexMisIDExtCut = new TCanvas("cvertexMisIDExtCut","");
  gStyle->SetOptStat(0.);
  cvertexMisIDExtCut->SetGridx(true);
  cvertexMisIDExtCut->SetGridy(true);
  hvertexMisIDExtCut->Divide(hvertexMisIDExtCut_total);
  hvertexMisIDExtCut->SetLineWidth(2);
  hvertexMisIDExtCut->SetLineColor(kBlue);
  hvertexMisIDExtCut->GetXaxis()->SetRangeUser(0,2000); 
  hvertexMisIDExtCut->GetYaxis()->SetRangeUser(0,1.); 
  hvertexMisIDExtCut->Draw();
  hvertexMisIDExtCut->SaveAs(Form("root_output/hvertexMisIDExtCut_%s.root",file->GetName()));
  cvertexMisIDExtCut->SaveAs("plots/vertexMisIDExtCut.eps");

  cout<<"End of this file"<<endl;
  
  if(!monoEnergy){
    delete hvertexXgood;
    delete hvertexXnhits;
    delete hvertexXdwall;
    delete hvertexXtowall;
    delete hvertexLongitudinalXTime;delete hvertexLongitudinal;delete hvertexOrtho;
    delete hvertex;delete hdirection;delete hnhits;
    delete hdwallXdwallrec;
    delete gdwall;
    delete gtowall;
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
  bool monoNRJ=true;
  bool energyOnly=true;
  double PC = 20;
  double DR = 1.0;//4.2;//8.4;

  char filename[256];
  sprintf(filename,"LEAF_%.0fPC_%.1fkHz",PC,DR);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  TApplication theApp("theApp", &argc, argv);
  
  double vertex_resolution, efficiency, migrationProbability;  
  
  if(monoNRJ){

    vertexResolution(10,vertex_resolution,efficiency,migrationProbability,monoNRJ,PC,DR);
    cout<<"Efficiency ="<<efficiency<<endl;
    cout<<"vertex Resolution ="<<vertex_resolution<<endl;

  }  

  else if(energyOnly){    

    const int nNRJ = 7;double event_energy[nNRJ] = {3,4,5,6,8,10,15};
    double vertex_resolution_energy[nNRJ];
    double efficiency_energy[nNRJ];
    double migrationProbability_energy[nNRJ];

    for(int ie=0;ie<nNRJ;ie++){
      vertexResolution(event_energy[ie],vertex_resolution,efficiency,migrationProbability,monoNRJ,PC,DR);
      cout<<"Efficiency ="<<efficiency<<", vertex Resolution ="<<vertex_resolution<<endl;
      vertex_resolution_energy[ie] = vertex_resolution;
      efficiency_energy[ie] = efficiency;
      migrationProbability_energy[ie] = migrationProbability;
    }
    
    TCanvas * cEnergy = new TCanvas("cEnergy","");
    TGraph * genergy;
    genergy = new TGraph(nNRJ,event_energy,vertex_resolution_energy);
    genergy->GetXaxis()->SetTitle("E_{#nu} (MeV)");
    genergy->GetYaxis()->SetTitle("Vertex resolution (cm)");
    genergy->GetXaxis()->SetRangeUser(3,15);
    genergy->GetYaxis()->SetRangeUser(0,300);
    genergy->SetMarkerStyle(20);
    genergy->SetLineColor(kBlue);
    genergy->SetMarkerColor(kBlue);
    genergy->Draw("ALP");
    genergy->SaveAs(Form("plots/genergy_%s.root",filename));
    cEnergy->SaveAs("plots/vertexEnergy.eps");
    
    TCanvas * cEfficiency = new TCanvas("cEfficiency","");
    TGraph * gefficiency;
    gefficiency = new TGraph(nNRJ,event_energy,efficiency_energy);
    gefficiency->GetXaxis()->SetTitle("E_{#nu} (MeV)");
    gefficiency->GetYaxis()->SetTitle("Signal efficiency");
    gefficiency->GetXaxis()->SetRangeUser(3,15);
    gefficiency->GetYaxis()->SetRangeUser(0,1.05);
    gefficiency->SetMarkerStyle(20);
    gefficiency->SetLineColor(kBlue);
    gefficiency->SetMarkerColor(kBlue);
    gefficiency->Draw("ALP");
    gefficiency->SaveAs(Form("plots/gefficiency_%s.root",filename));
    cEfficiency->SaveAs("plots/vertexEfficiency.eps");    


    TCanvas * cMigrationProbability_mono = new TCanvas("cMigrationProbability_mono","");
    TGraph * gmigrationProbability;
    gmigrationProbability = new TGraph(nNRJ,event_energy,migrationProbability_energy);
    gmigrationProbability->GetXaxis()->SetTitle("E_{#nu} (MeV)");
    gmigrationProbability->GetYaxis()->SetTitle("Signal migrationProbability");
    gmigrationProbability->GetXaxis()->SetRangeUser(3,15);
    gmigrationProbability->GetYaxis()->SetRangeUser(0,1.05);
    gmigrationProbability->SetMarkerStyle(20);
    gmigrationProbability->SetLineColor(kBlue);
    gmigrationProbability->SetMarkerColor(kBlue);
    gmigrationProbability->Draw("ALP");
    gmigrationProbability->SaveAs(Form("plots/gmigrationProbability_%s.root",filename));
    cMigrationProbability_mono->SaveAs("plots/vertexMigrationProbability.eps");    
  }

  else{
    const int nNRJ = 7;double event_energy[nNRJ] = {3,4,5,6,8,10,15};
    const int nPC = 3;double event_pc[nPC] = {20.,30.,40.};
    const int nDR = 4;double event_dr[nDR] = {4.2,6.0,8.4,10.};
    int color[nPC] = {600,416+2,1};//,632,416+2};
    int style[nDR] = {10,7,2,1};

    double vertex_resolution_energy[nPC][nDR][nNRJ];
    double signal_to_noise[nPC][nDR][nNRJ];
    double efficiency_energy[nPC][nDR][nNRJ];
    double migrationProbability_energy[nPC][nDR][nNRJ];
    dataID.clear();
    int counter=0;
    bool minimize=false;
    

    if(minimize){
            
      for(int ipc=0;ipc<nPC;ipc++){   
	for(int idr=0;idr<nDR;idr++){   
	  for(int ie=0;ie<nNRJ;ie++){
	    vertexResolution(event_energy[ie],vertex_resolution,efficiency,monoNRJ,event_pc[ipc],event_dr[idr]);
	    //cout<<"Efficiency ="<<efficiency<<", vertex Resolution ="<<vertex_resolution<<endl;
	    vertex_resolution_energy[ipc][idr][ie] = vertex_resolution;
	    efficiency_energy[ipc][idr][ie] = efficiency;
	    array <double,ntuple> data={event_pc[ipc],event_dr[ipc],event_energy[ie],vertex_resolution};
	    dataID.push_back(data);
	    cout<<"array <double,ntuple> data"<<counter<<"={"<<event_pc[ipc]<<","<<event_dr[idr]<<","<<event_energy[ie]<<","<<vertex_resolution<<"};dataID.push_back(data"<<counter<<");"<<endl;
	    counter++;
	  }
	}
      }
      double n[dimNoise]={0.};
      double s[dimSignal]={0.};
      minimizeNoise(n,s);
      cout<<"For points: noise: ";
      for(int i=0;i<dimNoise;i++) cout<<n[i]<<", ";
      cout<<" signal: ";
      for(int i=0;i<dimSignal;i++) cout<<s[i]<<", ";
      cout<<endl;
      
      calculateRMS(n,s,true);
    }
    else{
      double n[dimNoise]={970.062, 2.47523e-08, 1.77921};
      double s[dimSignal]={0.,0.,1.};
      for(int ipc=0;ipc<nPC;ipc++){   
	for(int idr=0;idr<nDR;idr++){   
	  for(int ie=0;ie<nNRJ;ie++){
	    vertexResolution(event_energy[ie],vertex_resolution,efficiency,migrationProbability,monoNRJ,event_pc[ipc],event_dr[idr]);
	    //cout<<"PMT type ="<<pmtType<<", Efficiency ="<<efficiency<<", vertex Resolution ="<<vertex_resolution<<endl;
	    vertex_resolution_energy[ipc][idr][ie] = vertex_resolution;
	    efficiency_energy[ipc][idr][ie] = efficiency;
	    migrationProbability_energy[ipc][idr][ie] = migrationProbability;
	    signal_to_noise[ipc][idr][ie]=SignalToNoise(event_pc[ipc],event_dr[idr],event_energy[ie],n,s);
	    array <double,ntuple> data={event_pc[ipc],event_dr[ipc],event_energy[ie],vertex_resolution};
	    dataID.push_back(data);
	    counter++;
	  }
	}
      }
    }
    
    TGraph * genergy[nPC][nDR];
    TGraph * gefficiency[nPC][nDR];
    TGraph * gmigrationProbability[nPC][nDR];
    counter=0;
    
    TCanvas * cEnergy = new TCanvas("cEnergy","");
    TCanvas * cEnergy_perDR[nDR];
    TCanvas * cEnergy_perPC[nPC];
    for(int ipc=0;ipc<nPC;ipc++){
      for(int idr=0;idr<nDR;idr++){
	if(ipc==0) cEnergy_perDR[idr] = new TCanvas(Form("cEnergy_perDR_%d",idr),"");
	if(idr==0) cEnergy_perPC[ipc] = new TCanvas(Form("cEnergy_perPC_%d",ipc),"");
	genergy[ipc][idr] = new TGraph(nNRJ,event_energy,vertex_resolution_energy[ipc][idr]);
	genergy[ipc][idr]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	genergy[ipc][idr]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	genergy[ipc][idr]->GetXaxis()->SetRangeUser(3,15);
	genergy[ipc][idr]->GetYaxis()->SetRangeUser(0,300);
	genergy[ipc][idr]->SetMarkerStyle(20);
	genergy[ipc][idr]->SetLineColor(color[ipc]);
	genergy[ipc][idr]->SetMarkerColor(color[ipc]);
	genergy[ipc][idr]->SetLineStyle(style[idr]);
	cEnergy->cd();
	if(counter==0) genergy[ipc][idr]->Draw("ALP");
	else genergy[ipc][idr]->Draw("sameLP");
	cEnergy_perDR[idr]->cd();
	if(ipc==0) genergy[ipc][idr]->Draw("ALP");
	else genergy[ipc][idr]->Draw("sameLP");
	cEnergy_perPC[ipc]->cd();
	if(idr==0) genergy[ipc][idr]->Draw("ALP");
	else genergy[ipc][idr]->Draw("sameLP");
	//genergy[ipc][idr]->SaveAs(Form("plots/genergy[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	counter++;
      }
    }
    
    cEnergy->SaveAs("plots/vertexEnergy.eps");
    for(int ipc=0;ipc<nPC;ipc++) cEnergy_perPC[ipc]->SaveAs(Form("plots/vertexEnergy_PC%d.eps",ipc));
    for(int idr=0;idr<nDR;idr++) cEnergy_perDR[idr]->SaveAs(Form("plots/vertexEnergy_DR%d.eps",idr));

    TCanvas * cSignalToNoise = new TCanvas("cSignalToNoise","");
    TGraph * gsignal_to_noise[nPC][nDR];
    counter=0;

    for(int ipc=0;ipc<nPC;ipc++){   
      for(int idr=0;idr<nDR;idr++){
	gsignal_to_noise[ipc][idr] = new TGraph(nNRJ,signal_to_noise[ipc][idr],vertex_resolution_energy[ipc][idr]);
	gsignal_to_noise[ipc][idr]->GetXaxis()->SetTitle("S/#sqrt(N)");
	gsignal_to_noise[ipc][idr]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	//gsignal_to_noise[ipc][idr]->GetXaxis()->SetRangeUser(3,15);
	gsignal_to_noise[ipc][idr]->GetYaxis()->SetRangeUser(0,300);
	gsignal_to_noise[ipc][idr]->SetMarkerStyle(20);
	gsignal_to_noise[ipc][idr]->SetLineColor(color[ipc]);
	gsignal_to_noise[ipc][idr]->SetMarkerColor(color[ipc]);
	gsignal_to_noise[ipc][idr]->SetLineStyle(style[idr]);
	if(counter==0) gsignal_to_noise[ipc][idr]->Draw("ALP");
	else gsignal_to_noise[ipc][idr]->Draw("sameLP");
	//gsignal_to_noise[ipc][idr]->SaveAs(Form("plots/gsignal_to_noise[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	counter++;
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
	gefficiency[ipc][idr] = new TGraph(nNRJ,event_energy,efficiency_energy[ipc][idr]);
	gefficiency[ipc][idr]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	gefficiency[ipc][idr]->GetYaxis()->SetTitle("Vertex resolution (cm)");
	gefficiency[ipc][idr]->GetXaxis()->SetRangeUser(3,15);
	gefficiency[ipc][idr]->GetYaxis()->SetRangeUser(0,1);
	gefficiency[ipc][idr]->SetMarkerStyle(20);
	gefficiency[ipc][idr]->SetLineColor(color[ipc]);
	gefficiency[ipc][idr]->SetMarkerColor(color[ipc]);
	gefficiency[ipc][idr]->SetLineStyle(style[idr]);
	cEfficiency->cd();
	if(counter==0) gefficiency[ipc][idr]->Draw("ALP");
	else gefficiency[ipc][idr]->Draw("sameLP");
	cEfficiency_perDR[idr]->cd();
	if(ipc==0) gefficiency[ipc][idr]->Draw("ALP");
	else gefficiency[ipc][idr]->Draw("sameLP");
	cEfficiency_perPC[ipc]->cd();
	if(idr==0) gefficiency[ipc][idr]->Draw("ALP");
	else gefficiency[ipc][idr]->Draw("sameLP");
	//gefficiency[ipc][idr]->SaveAs(Form("plots/gefficiency[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	counter++;
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
	gmigrationProbability[ipc][idr] = new TGraph(nNRJ,event_energy,migrationProbability_energy[ipc][idr]);
	gmigrationProbability[ipc][idr]->GetXaxis()->SetTitle("E_{#nu} (MeV)");
	gmigrationProbability[ipc][idr]->GetYaxis()->SetTitle("Migration probability");
	gmigrationProbability[ipc][idr]->GetXaxis()->SetRangeUser(3,15);
	gmigrationProbability[ipc][idr]->GetYaxis()->SetRangeUser(0,0.4);
	gmigrationProbability[ipc][idr]->SetMarkerStyle(20);
	gmigrationProbability[ipc][idr]->SetLineColor(color[ipc]);
	gmigrationProbability[ipc][idr]->SetMarkerColor(color[ipc]);
	gmigrationProbability[ipc][idr]->SetLineStyle(style[idr]);
	cMigrationProbability->cd();
	if(counter==0) gmigrationProbability[ipc][idr]->Draw("ALP");
	else gmigrationProbability[ipc][idr]->Draw("sameLP");
	cMigrationProbability_perDR[idr]->cd();
	if(ipc==0) gmigrationProbability[ipc][idr]->Draw("ALP");
	else gmigrationProbability[ipc][idr]->Draw("sameLP");
	cMigrationProbability_perPC[ipc]->cd();
	if(idr==0) gmigrationProbability[ipc][idr]->Draw("ALP");
	else gmigrationProbability[ipc][idr]->Draw("sameLP");
	//gmigrationProbability[ipc][idr]->SaveAs(Form("plots/gmigrationProbability[ipc][idr]_pmtType%d_%s.root",pmtType,filename));
	counter++;
      }
    }
    cMigrationProbability->SaveAs("plots/vertexMigrationProbability.eps");
    for(int ipc=0;ipc<nPC;ipc++) cMigrationProbability_perPC[ipc]->SaveAs(Form("plots/vertexMigrationProbability_PC%d.eps",ipc));
    for(int idr=0;idr<nDR;idr++) cMigrationProbability_perDR[idr]->SaveAs(Form("plots/vertexMigrationProbability_DR%d.eps",idr));
  }
  
  theApp.Run();
  return 0;
}

  
