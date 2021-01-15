void ProducePDF(){
  TString ofname = "timePDF.root";
  TString ifname = "Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_10MeV.root"; // Originally used in B.Quilain codes

  int opt;

  const int nFiles=1;
  TFile * _f[nFiles];
  const int nPMTtypes = 2;
  double DRTotalPerNS[nPMTtypes];

  //for(int i=0;i<nPMTtypes;i++){
  //if(i==0) DRTotalPerNS[i]=8400*2e4*1e-9;
  //else if (i==1) DRTotalPerNS[i]=100*19*5e3*1e-9;
  //}
  bool isNew=true;//false;//true;//false;
  bool isNew2=false;
  
  int colorScale[nPMTtypes];
  for(int i=0;i<nPMTtypes;i++){
    if(i==0) colorScale[i] = 4;
    else if(i==1) colorScale[i] = 2;
    else colorScale[i] = 3;
  }

  int colorScale2[4];
  for(int i=0;i<4;i++){
    if(i==0) colorScale2[i] = 2;
    else if(i==1) colorScale2[i] = 4;
    else if(i==2) colorScale2[i] = 416+2;
    else colorScale2[i] = 1;
  }

  TLegend * l = new TLegend(0.65,0.65,0.89,0.89);
  l->SetLineColor(0);

  TFile * fOut = new TFile(ofname,"recreate");
  
  TH1D * ChargeProfile[nFiles][nPMTtypes];TCanvas * cChargeProfile;TH1D * ratioChargeProfile[nFiles][nPMTtypes];
  TH1D * ChargePerPMT[nFiles][nPMTtypes];TCanvas *cChargePerPMT;
  TH1D * TotalCharge[nFiles][nPMTtypes];TCanvas *cTotalCharge; TH1D * TotalHit[nFiles][nPMTtypes]; TCanvas *cTotalHit;

  TH1D * TimeProfile[nFiles][nPMTtypes];TCanvas * cTimeProfile;
  TH1D * TimeTOFProfile[nFiles][nPMTtypes];TCanvas * cTimeTOFProfile;
  TH1D * HitTimeTOFProfile[nFiles][nPMTtypes];TCanvas * cHitTimeTOFProfile;
  TH1D * HitTimeTOFDR[nFiles][nPMTtypes];
  
  TH2D * TimeTOFProfileXTOF[nFiles][nPMTtypes];TCanvas * cTimeTOFProfileXTOF;

  const int nGroupsPMTs = 3;
  TH3D * angularResponseXPMT_hits[nFiles][nPMTtypes][nGroupsPMTs];TH3D * clone_angularResponseXPMT_hits[nFiles][nPMTtypes][nGroupsPMTs];
  TH3D * angularResponseXPMT_allPMTs[nFiles][nPMTtypes][nGroupsPMTs];TH3D * clone_angularResponseXPMT_allPMTs[nFiles][nPMTtypes][nGroupsPMTs];
  const int nBinsDistance=1;  
  const int nBinsTheta=1;//8
  //const int nBinsPhi=1;
  bool killLowStatBin=true;//Remove statistical fluctuations to generate angular PDF
  
  TCanvas * cAngularResponsePMT_hits[nPMTtypes][nGroupsPMTs][nBinsDistance];
  TH2D *angular2DResponse_hits[nPMTtypes][nGroupsPMTs][nBinsDistance]; TH2D *angular2DResponse_allPMTs[nPMTtypes][nGroupsPMTs][nBinsDistance];
  TGraph2D * gAngular2DResponse[nPMTtypes][nGroupsPMTs][nBinsDistance];
  
  TCanvas * cAngularXDistResponsePMT_hits[nPMTtypes][nGroupsPMTs][nBinsDistance];
  TH2D *angularXDist2DResponse_hits[nPMTtypes][nGroupsPMTs][nBinsDistance]; TH2D *angularXDist2DResponse_allPMTs[nPMTtypes][nGroupsPMTs][nBinsDistance];

  TCanvas * cDistResponsePMT_hits[nPMTtypes][nGroupsPMTs];
  TH1D * hDistResponsePMT_hits[nPMTtypes][nGroupsPMTs][nBinsTheta]; TH1D * hDistResponsePMT_allPMTs[nPMTtypes][nGroupsPMTs][nBinsTheta];
  TF1 * fDistResponsePMT[nPMTtypes][nGroupsPMTs][nBinsTheta];
  TGraph * gDistResponsePMT[nPMTtypes][nGroupsPMTs][nBinsTheta];
  TSpline3 * sDistResponsePMT[nPMTtypes][nGroupsPMTs][nBinsTheta];
  
  //TCanvas * cAngularResponsePMT_hits_1D[nPMTtypes][nBinsDistance];
 
  TH1D *hProjection_1D;TH1D *hProjection2_1D;
  TH1D *angular2DResponse_hits_1D[nPMTtypes][nGroupsPMTs];

  TH3D * ChargeProfile2DXTheta_hits[nFiles][nPMTtypes];TCanvas * cChargeProfile2DXTheta_hits[nPMTtypes][nGroupsPMTs];
  TH3D * ChargeProfile2DXTheta_allPMTs[nFiles][nPMTtypes];
  TH2D * hChargeProjection[nPMTtypes][nGroupsPMTs];
  TH2D * hChargeProjection2[nPMTtypes][nGroupsPMTs];
  TH2D * ChargeProfile2D_hits[nPMTtypes][nGroupsPMTs];

  //TH3D * ChargeProfile2DXPMT[nFiles][nPMTtypes];TCanvas * cChargeProfile2DXPMT;

  TH2D * ChargeProfileXdWall[nFiles][nPMTtypes];TCanvas * cChargeProfileXdWall;//TH1D * ratioChargeProfile[nFiles][nPMTtypes];

  const int nBinsTOF = 10;
  TH1D * TimeTOFProfileXTOF_1D[nFiles][nPMTtypes][nBinsTOF];TCanvas * cTimeTOFProfileXTOF_1D;
 
  TF1 * gausExpoConv[nFiles][nPMTtypes];
  TF1 * ExpoQueue[nFiles][nPMTtypes];
  TF1 * gausExpoConv_slice[nFiles][nPMTtypes][nBinsTOF];
  TSpline3 * splineExpoConv[nFiles][nPMTtypes];
  TSpline3 * splineExpoQueue[nFiles][nPMTtypes];
  TSpline3 * splineDR[nFiles][nPMTtypes];
  double limitFitGausExpo[nPMTtypes];
  for(int i=0;i<nPMTtypes;i++){
    if(i==0) limitFitGausExpo[i] = 800;//100
    else limitFitGausExpo[i]=800;//100
  }
  TGraph * graphExpoQueue[nFiles][nPMTtypes];
  TGraph * graphDR[nFiles][nPMTtypes];
  TH3D * hPMTDirectionality[nPMTtypes];

  
  for(int f=0;f<nFiles;f++){

    //_f[f] = new TFile("test6_ampute.root","read");
    //_f[f] = new TFile("test6_uniform.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT_10MeV.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374100Hz_3MeV_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_noDN_digitized_v2.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("test.root","read");
    //_f[f] = new TFile("pdfHE_2.root","read");
    //_f[f] = new TFile("pdfHE_500MeVelectron.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized_fullStat_50Hz.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_3MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_8MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_3MeV_digitized_fullStat_noDR.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_2.2MeVgamma_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_10mPMT14374_2.2MeVgamma_digitized.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_10mPMT14374_15MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_10mPMT14374_10MeV_digitized_fullStat.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_digitized_fullStat_all.root","read");

    //_f[f] = new TFile("Hybrid_20BAL_5mPMT14374_10MeV_noDR_digitized_fullStat_all.root","read");
    //if(f==0) _f[f] = new TFile("test_hk20pc_4.2kHzbl_e10.root","read");
    //if(f==0) _f[f] = new TFile("Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_10MeV.root","read");
    //if(f==0) _f[f] = new TFile("Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_10MeV_nodr.root","read");
    //if(f==0) _f[f] = new TFile("Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_3MeV.root","read");
    //if(f==0) _f[f] = new TFile("10MeV_nodr.root","read");
    //else _f[f] = new TFile("test_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e10.root","read");
    //_f[f] = new TFile("gamma2.2MeV_5kmPMT.root","read");
    //_f[f] = new TFile("gamma2.2MeV_10kmPMT.root","read");
    //_f[f] = new TFile("inputs_10MeV_merged.root","read");
    //_f[f] = new TFile("Hybrid_20BAL_5mPMT_digitized.root","read");
    //_f[f] = new TFile("test7_uniform.root","read");
    //_f[f] = new TFile("test6_ampute_uniform.root","read");
    //_f[f] = new TFile("test6_ampute_digitized.root","read");

    if(f==0) _f[f] = new TFile(ifname,"read");

    double nEvents=0;
    double events[nPMTtypes];
    
    for(int i=0;i<nPMTtypes;i++){
      
      ChargePerPMT[f][i] = (TH1D*) _f[f]->Get(Form("ChargePerPMT_pmtType%d",i));
      ChargePerPMT[f][i]->SetLineWidth(2);
      ChargePerPMT[f][i]->SetLineColor(colorScale[i]);
      ChargePerPMT[f][i]->GetXaxis()->SetTitle("Q (p.e)");
      ChargePerPMT[f][i]->GetYaxis()->SetTitle("Number of events (normalized)");
      if(f==0){
	if(i==0) l->AddEntry(ChargePerPMT[f][i],"B&L");
	else if(i==1) l->AddEntry(ChargePerPMT[f][i],"multi-PMT");
      }

      ChargeProfile[f][i] = (TH1D*) _f[f]->Get(Form("ChargeProfile_pmtType%d",i));
      ChargeProfile[f][i]->SetLineWidth(2);
      ChargeProfile[f][i]->SetLineColor(colorScale[i]);
      ChargeProfile[f][i]->GetXaxis()->SetTitle("Angle with particle direction (#circ)");
      ChargeProfile[f][i]->GetYaxis()->SetTitle("Charge (p.e)");

      if(isNew2){
	ChargeProfileXdWall[f][i] = (TH2D*) _f[f]->Get(Form("ChargeProfileXdWall_pmtType%d",i));
	ChargeProfileXdWall[f][i]->GetXaxis()->SetTitle("dwall (cm)");
	ChargeProfileXdWall[f][i]->GetYaxis()->SetTitle("Angle with particle direction (#circ)");
      }
      
      TotalCharge[f][i] = (TH1D*) _f[f]->Get(Form("TotalCharge_pmtType%d",i));
      TotalCharge[f][i]->SetLineWidth(2);
      TotalCharge[f][i]->SetLineColor(colorScale[i]);
      TotalCharge[f][i]->GetXaxis()->SetTitle("Q (p.e)");
      TotalCharge[f][i]->GetYaxis()->SetTitle("Number of events (normalized)");
      nEvents += TotalCharge[f][i]->Integral();

      events[i] = TotalCharge[f][i]->Integral();
	
      TotalHit[f][i] = (TH1D*) _f[f]->Get(Form("TotalHit_pmtType%d",i));
      TotalHit[f][i]->SetLineWidth(2);
      TotalHit[f][i]->SetLineColor(colorScale[i]);
      TotalHit[f][i]->GetXaxis()->SetTitle("Q (p.e)");
      TotalHit[f][i]->GetYaxis()->SetTitle("Number of events (normalized)");

      TimeProfile[f][i] = (TH1D*) _f[f]->Get(Form("TimeProfile_pmtType%d",i));
      TimeProfile[f][i]->SetLineWidth(2);
      TimeProfile[f][i]->SetLineColor(colorScale[i]);
      TimeProfile[f][i]->GetXaxis()->SetTitle("Time (ns)");
      TimeProfile[f][i]->GetYaxis()->SetTitle("Charge (p.e) (normalized)");

      TimeTOFProfile[f][i] = (TH1D*) _f[f]->Get(Form("TimeTOFProfile_pmtType%d",i));
      TimeTOFProfile[f][i]->SetLineWidth(2);
      TimeTOFProfile[f][i]->SetLineColor(colorScale[i]);
      TimeTOFProfile[f][i]->GetXaxis()->SetTitle("Time - TOF (ns)");
      TimeTOFProfile[f][i]->GetYaxis()->SetTitle("Charge (p.e) (normalized)");

      HitTimeTOFProfile[f][i] = (TH1D*) _f[f]->Get(Form("HitTimeTOFProfile_pmtType%d",i));
      HitTimeTOFProfile[f][i]->SetLineWidth(2);
      HitTimeTOFProfile[f][i]->SetLineColor(colorScale[i]);
      HitTimeTOFProfile[f][i]->GetXaxis()->SetTitle("HitTime - TOF (ns)");
      HitTimeTOFProfile[f][i]->GetYaxis()->SetTitle("nhits (normalized)");

      HitTimeTOFDR[f][i] = (TH1D*)  HitTimeTOFProfile[f][i]->Clone(Form("HitTimeTOFDR_pmtType%d",i));
      HitTimeTOFDR[f][i]->Reset();

      
      TimeTOFProfileXTOF[f][i] = (TH2D*) _f[f]->Get(Form("TimeTOFProfileXTOF_pmtType%d",i));
      //TimeTOFProfileXTOF[f][i]->SetLineWidth(2);
      TimeTOFProfileXTOF[f][i]->GetXaxis()->SetTitle("Time - TOF (ns)");
      TimeTOFProfileXTOF[f][i]->GetYaxis()->SetTitle("TOF (ns)");
      
      for(int p=0;p<nGroupsPMTs;p++){
	angularResponseXPMT_hits[f][i][p] = (TH3D*) _f[f]->Get(Form("angularResponseXPMT_hits_pmtType%d_pmtgroup%d",i,p));
	angularResponseXPMT_allPMTs[f][i][p] = (TH3D*) _f[f]->Get(Form("angularResponseXPMT_allPMTs_pmtType%d_pmtgroup%d",i,p));
	angularResponseXPMT_hits[f][i][p]->Sumw2();
	angularResponseXPMT_allPMTs[f][i][p]->Sumw2();
	clone_angularResponseXPMT_hits[f][i][p] = (TH3D*) angularResponseXPMT_hits[f][i][p]->Clone(Form("angularResponseXPMT_hits_pmtType%d_pmtgroup%d",i,p));
	clone_angularResponseXPMT_allPMTs[f][i][p] = (TH3D*) angularResponseXPMT_allPMTs[f][i][p]->Clone(Form("angularResponseXPMT_allPMTs_pmtType%d_pmtgroup%d",i,p));
	clone_angularResponseXPMT_hits[f][i][p]->Sumw2();
	clone_angularResponseXPMT_allPMTs[f][i][p]->Sumw2();

	//angularResponseXPMT_hits[f][i][p]->Rebin3D(1,3,1);
	//angularResponseXPMT_allPMTs[f][i][p]->Rebin3D(1,3,1);
	angularResponseXPMT_hits[f][i][p]->Rebin3D(1,18,1);
	angularResponseXPMT_allPMTs[f][i][p]->Rebin3D(1,18,1);
	//angularResponseXPMT_hits[f][i][p]->Rebin3D(1,2,1);
	//angularResponseXPMT_allPMTs[f][i][p]->Rebin3D(1,2,1);
	clone_angularResponseXPMT_hits[f][i][p]->Rebin3D(2,18,5);
	clone_angularResponseXPMT_allPMTs[f][i][p]->Rebin3D(2,18,5);
      }
      ChargeProfile2DXTheta_hits[f][i] = (TH3D*) _f[f]->Get(Form("ChargeProfile2DXTheta_hits_pmtType%d",i));
      ChargeProfile2DXTheta_allPMTs[f][i] = (TH3D*) _f[f]->Get(Form("ChargeProfile2DXTheta_allPMTs_pmtType%d",i));
      cout<<"Ola"<<endl;
    }
    

    fOut->cd();	      
    TH2D * hProjection;
    TH2D * hProjection2;
    for(int i=0;i<nPMTtypes;i++){
      for(int p=0;p<nGroupsPMTs;p++){
	
	//Prepare the distance-only response:
	if(f==0) cDistResponsePMT_hits[i][p] = new TCanvas(Form("cDistResponsePMT_hits_pmtType%d_group%d",i,p),Form("cDistResponsePMT_hits_pmtType%d_group%d",i,p));
	
	//clone_angularResponseXPMT_hits[f][i][p]->Rebin3D(2,18,1);
	//clone_angularResponseXPMT_allPMTs[f][i][p]->Rebin3D(2,18,1);
	
	for(int ibinx=0;ibinx<nBinsTheta;ibinx++){
	  //for(int ibiny=0;ibiny<nBinsPhi;ibiny++){
	    cout<<"Theta = "<<ibinx<<endl;
	    //cout<<"Theta = "<<ibinx<<", Phi = "<<ibiny<<endl;
	    
	    clone_angularResponseXPMT_hits[f][i][p]->GetZaxis()->SetRange(ibinx+1,ibinx+1);
	    clone_angularResponseXPMT_allPMTs[f][i][p]->GetZaxis()->SetRange(ibinx+1,ibinx+1);
	    //angularResponseXPMT_hits[f][i][p]->GetXaxis()->SetRange(ibinx+1,ibinx+1);
	    //angularResponseXPMT_allPMTs[f][i][p]->GetXaxis()->SetRange(ibinx+1,ibinx+1);
	    
	    hProjection_1D = (TH1D*) clone_angularResponseXPMT_hits[f][i][p]->Project3D("x");
	    hProjection2_1D = (TH1D*) clone_angularResponseXPMT_allPMTs[f][i][p]->Project3D("x");
	    //hProjection->Rebin2D(2,2);
	    //hProjection2->Rebin2D(2,2);
	    hProjection_1D->Divide(hProjection2_1D);
	    
	    hDistResponsePMT_hits[i][p][ibinx] = (TH1D*) hProjection_1D->Clone(Form("hDistResponsePMT_hits_pmtType%d_pmt%d_theta%d",i,p,ibinx));
	    for(int ib=0;ib<hDistResponsePMT_hits[i][p][ibinx]->GetNbinsX();ib++){
	      double content=hDistResponsePMT_hits[i][p][ibinx]->GetBinContent(ib+1);
	      double error=hDistResponsePMT_hits[i][p][ibinx]->GetBinError(ib+1);
	      double dist=hDistResponsePMT_hits[i][p][ibinx]->GetBinCenter(ib+1);
	      //double value=
	      //hDistResponsePMT_hits[i][p][ibinx]->SetBinContent(ib+1,content*pow(dist,4));
	      //hDistResponsePMT_hits[i][p][ibinx]->SetBinError(ib+1,error*pow(dist,4));
	    }
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"x > [6] ? [0] + [1]*x + [2]*(x*x) + [3]*(x*x*x) + [4]*(x*x*x*x) : [7] + 8*x",0,1e5);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"expo",0,50);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"pol1",0,50);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"pol4",0,1e5);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"pol2",0,50);
	    double cutValue=100;
	    //TF1 * ftemp = new TF1("temp","pol2",0,50);
	    //hDistResponsePMT_hits[i][p][ibinx]->Fit("temp","R");
	    //double parFirstFit[3]=ftemp->GetParameters();
	    fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"x > [6] ? [0] + [1]/x + [2]/(x*x) + [3]/(x*x*x) + [4]/(x*x*x*x) : [7] + [8]*x + [9]*x*x",0,1e5);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"x > [6] ? [0] + [1]/x + [2]/(x*x) + [3]/(x*x*x) + [4]/(x*x*x*x) : [7] + x*[8]",0,1e5);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx),"x > [6] ? [0] + [1]/x + [2]/(x*x) + [3]/(x*x*x) + [4]/(x*x*x*x) : [7]*TMath::Exp(-x/[8])",0,1e5);
	    //fDistResponsePMT[i][p][ibinx] = new TF1(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d_phi%d",i,p,ibinx,ibiny),"[0] + [1]/x + [2]/(x*x) + [3]/(x*x*x) + [4]/(x*x*x*x) + [5]/(x*x*x*x*x) + [6]/(x*x*x*x*x*x)",0,1e5);
	    fDistResponsePMT[i][p][ibinx]->FixParameter(6,cutValue);//Around 30 cm
	    //for(int ip=0;ip<3;ip++){
	    //fDistResponsePMT[i][p][ibinx]->FixParameter(7+ip,ftemp->GetParameter(ip));//Around 30 cm
	    //}
	    //delete ftemp;
	    
	    //fDistResponsePMT[i][p][ibinx]->SetParameter(7,hDistResponsePMT_hits[i][p][ibinx]->GetBinContent(hDistResponsePMT_hits[i][p][ibinx]->FindBin(40)));//Around 30 cm
	    /*
	    fDistResponsePMT[i][p][ibinx]->SetParameter(0,3.011e5);
	    fDistResponsePMT[i][p][ibinx]->SetParameter(1,-4.417e4);
	    fDistResponsePMT[i][p][ibinx]->SetParameter(2,1.463e3);
	    fDistResponsePMT[i][p][ibinx]->SetParameter(3,-8.001e-2);
	    fDistResponsePMT[i][p][ibinx]->SetParameter(4,1.415e-4);
	    */	    
	    fDistResponsePMT[i][p][ibinx]->SetLineColor(2);
	    fDistResponsePMT[i][p][ibinx]->SetLineWidth(2);
	    hDistResponsePMT_hits[i][p][ibinx]->Fit(Form(Form("fDistResponsePMT_pmtType%d_pmt%d_theta%d",i,p,ibinx)));
	    //hDistResponsePMT_allPMTs[i][p][ibinx] = (TH1D*) hProjection2_1D->Clone(Form("hDistResponsePMT_allPMTs_pmtType%d_pmt%d_theta%d_phi%d",i,p,ibinx,ibiny));
	    cDistResponsePMT_hits[i][p]->cd(1);
	    hDistResponsePMT_hits[i][p][ibinx]->SetLineColor(ibinx+1);
	    hDistResponsePMT_hits[i][p][ibinx]->SetLineWidth(2);
	    if(ibinx==0) hDistResponsePMT_hits[i][p][ibinx]->Draw("E1hist");
	    else hDistResponsePMT_hits[i][p][ibinx]->Draw("E1histsame");
	    fDistResponsePMT[i][p][ibinx]->Draw("same");
	    fDistResponsePMT[i][0][ibinx]->Draw("same");
	    if(p==0) fDistResponsePMT[i][p][ibinx]->Write(Form(Form("fDistResponsePMT_pmtType%d",i)));
	    //hProjection[i][p]->Draw("colz");	  
	    //angularXDist2DResponse_hits[i][p][ibinx]->Draw("colz");	  
	    //cAngularXDistResponsePMT_hits[i][p][ibinx]->cd(2);
	    //angularXDist2DResponse_allPMTs[i][p]->Draw("colz");	  
	    //hProjection2[i][p]->Draw("colz");	  
	    //cAngularXDistResponsePMT_hits[i][p]->SetLogz();
	    //cAngularXDistResponsePMT_hits[i][p][ibinx]->SaveAs(Form("plots/cAngularXDistResponsePMT2D_pmtType%d_group%d_distance%d.eps",i,p,ibinx));
	    //angularXDist2DResponse_hits[i][p]->Write(Form("hPMTDirectionality_%d_%d_%d",f,i,p));
	    //}
	}



	
	//Second, we will integrate out the dependency on distance by assuming it goes as PMT-vertex distance squared.
	for(int ibiny=0;ibiny<angularResponseXPMT_hits[f][i][p]->GetNbinsY();ibiny++){
	  for(int ibinz=0;ibinz<angularResponseXPMT_hits[f][i][p]->GetNbinsZ();ibinz++){
	    for(int ibinx=0;ibinx<angularResponseXPMT_hits[f][i][p]->GetNbinsX();ibinx++){
	      double distance = angularResponseXPMT_hits[f][i][p]->GetXaxis()->GetBinCenter(ibinx+1);
	      double content_hits = angularResponseXPMT_hits[f][i][p]->GetBinContent(ibinx+1,ibiny+1,ibinz+1);
	      //double correction = pow(distance,2);
	      double correction = 1./fDistResponsePMT[i][p][0]->Eval(distance);
	      //if(i==1 && content_hits!=0) cout<<"Content = "<<content_hits<<", distance = "<<distance<<", value = "<<1./fDistResponsePMT[i][p][0]->Eval(distance)<<endl;
	      //content_hits*=1./fDistResponsePMT[i][p][0]->Eval(distance);
	      content_hits*=correction;
	      angularResponseXPMT_hits[f][i][p]->SetBinContent(ibinx+1,ibiny+1,ibinz+1,content_hits);
	      double error_hits = angularResponseXPMT_hits[f][i][p]->GetBinError(ibinx+1,ibiny+1,ibinz+1);
	      angularResponseXPMT_hits[f][i][p]->SetBinError(ibinx+1,ibiny+1,ibinz+1,error_hits*correction);
	      //double content_all = angularResponseXPMT_allPMTs[f][i][p]->GetBinContent(ibinx+1,ibiny+1,ibinz+1);
	      //content_all*=pow(distance,2);
	      //angularResponseXPMT_allPMTs[f][i][p]->SetBinContent(ibinx+1,ibiny+1,ibinz+1,content_all);
	    }
	  }
	}


	
	//Prepare the PMT angular theta X phi response
	for(int ibinx=0;ibinx<nBinsDistance;ibinx++){
	  //cout<<"distance bin = "<<ibinx<<endl;
	  if(f==0){
	    cAngularResponsePMT_hits[i][p][ibinx] = new TCanvas(Form("cAngularResponsePMT_hits_pmtType%d_group%d_distance%d",i,p,ibinx),Form("cAngularResponsePMT_hits_pmtType%d_group%d_distance%d",i,p,ibinx));
	    //cAngularResponsePMT_hits[i][p][ibinx]->Divide(2,1);
	  }
	
	  //int nbins = angularResponseXPMT_hits[f][i][p]->GetNbinsX();
	  //int binCutValue = angularResponseXPMT_hits[f][i][p]->GetXaxis()->FindBin(cutValue);
	  //angularResponseXPMT_hits[f][i][p]->GetXaxis()->SetRange(0,binCutValue);
	  //angularResponseXPMT_allPMTs[f][i][p]->GetXaxis()->SetRange(0,binCutValue);
	  //angularResponseXPMT_hits[f][i][p]->GetXaxis()->SetRange(ibinx+1,ibinx+1);
	  //angularResponseXPMT_allPMTs[f][i][p]->GetXaxis()->SetRange(ibinx+1,ibinx+1);
	  //angularResponseXPMT_hits[f][i][p]->GetXaxis()->SetRange(1,10);
	  //angularResponseXPMT_allPMTs[f][i][p]->GetXaxis()->SetRange(1,10);
	  hProjection = (TH2D*) angularResponseXPMT_hits[f][i][p]->Project3D("zy");
	  hProjection2 = (TH2D*) angularResponseXPMT_allPMTs[f][i][p]->Project3D("zy");
	  hProjection->Rebin2D(1,1);
	  hProjection2->Rebin2D(1,1);
	  hProjection->Divide(hProjection2);
	  
	  angular2DResponse_hits[i][p][ibinx] = (TH2D*) hProjection->Clone(Form("angular2DResponse_hits_pmtType%d_pmt%d_distance%d",i,p,ibinx));
	  angular2DResponse_allPMTs[i][p][ibinx] = (TH2D*) hProjection2->Clone(Form("angular2DResponse_allPMTs_pmtType%d_pmt%d_distance%d",i,p,ibinx));
	  cAngularResponsePMT_hits[i][p][ibinx]->cd(1);
	  //hProjection[i][p]->Draw("colz");

	  
	  //How to prepare the PDF:
	  //1. Replace the bin where the probability is 0 to a very small probability (0 which will just makes the likelihood too much discontinuous, especially, when adding DR)
	  //a. Look for the smallest bin content which is not 0
	  
	  //We will first kill bins with large stat. fluctuations. It is not done cleanly here, but the effect is the same. We will basically set their value to 0.
	  for(int ix=0;ix<angular2DResponse_hits[i][p][ibinx]->GetNbinsX();ix++){
	    for(int iy=0;iy<angular2DResponse_hits[i][p][ibinx]->GetNbinsY();iy++){
	      double content=angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	      double error=angular2DResponse_hits[i][p][ibinx]->GetBinError(ix+1,iy+1);
	      double rel_error=error;
	      if(content!=0) rel_error/=content;
	      //cout<<"Error of bin "<<ix<<", "<<iy<<" is "<<error<<" vs "<<content<<endl;
	      if(rel_error>1e-1 && killLowStatBin) angular2DResponse_hits[i][p][ibinx]->SetBinContent(ix+1,iy+1,0);
	    }
	  }

	  /*double * min = new double[angular2DResponse_hits[i][p][ibinx]->GetNbinsX()];
	  
	  for(int ix=0;ix<angular2DResponse_hits[i][p][ibinx]->GetNbinsX();ix++){
	  min[ix]=angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,1);
	    for(int iy=0;iy<angular2DResponse_hits[i][p][ibinx]->GetNbinsY();iy++){
	      double content=angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	      if(content!=0 && content<min[ix]) min[ix]=content;
	    }
	  }
	  */
	  
	  double min = 1e9;
	  for(int ix=0;ix<angular2DResponse_hits[i][p][ibinx]->GetNbinsX();ix++){
	    for(int iy=0;iy<angular2DResponse_hits[i][p][ibinx]->GetNbinsY();iy++){
	      double content=angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	      if(content!=0 && content<min) min=content;
	    }
	  }
	  
	  /*	
	  //For each phi, take the last non zero theta bin
	  //Amd evaluate the gradient of the decrese
	  for(int ix=0;ix<angular2DResponse_hits[i][p][ibinx]->GetNbinsX();ix++){
	    int iy = angular2DResponse_hits[i][p][ibinx]->GetNbinsY()-1;
	    double content = angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	    while(content==0 && iy>0){
	      iy--;
	      content = angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	      cout<<content<<endl;
	    }
	    double position = angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinCenter(iy+1);
	    double grad = (content-0.)/(position - 180.);
	    double init = content;
	    cout<<"Phi = "<<ix<<", theta = "<<position<<", content = "<<content<<endl;
	    
	    //Then, apply it.
	    for(int iy2=iy+1;iy2<angular2DResponse_hits[i][p][ibinx]->GetNbinsY();iy2++){
	      double pos =  angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinCenter(iy2+1);
	      double value = grad*(pos-position) + content;
	      angular2DResponse_hits[i][p][ibinx]->SetBinContent(ix+1,iy2+1,value);
	    }
	  }
	  */
		      
	  //b. Replace all bins = 0 by the smallest value
	  gAngular2DResponse[i][p][ibinx] = new TGraph2D();
	  int ig=0;
	  for(int ix=0;ix<angular2DResponse_hits[i][p][ibinx]->GetNbinsX();ix++){
	    for(int iy=0;iy<angular2DResponse_hits[i][p][ibinx]->GetNbinsY();iy++){
	      double content=angular2DResponse_hits[i][p][ibinx]->GetBinContent(ix+1,iy+1);
	      //if(content==0) content=min[ix];
	      if(content==0) content=min;
	      angular2DResponse_hits[i][p][ibinx]->SetBinContent(ix+1,iy+1,content);//=content;
	      //For points on the edges of the histograms, we will shift them to their more extreme edge to avoid TGraph2D to have e.g. no value for interpolation at (0,0)
	      double xposition = angular2DResponse_hits[i][p][ibinx]->GetXaxis()->GetBinCenter(ix+1);
	      if(ix == 0) xposition = angular2DResponse_hits[i][p][ibinx]->GetXaxis()->GetBinLowEdge(ix+1);
	      else if(ix == (angular2DResponse_hits[i][p][ibinx]->GetNbinsX()-1)) xposition = angular2DResponse_hits[i][p][ibinx]->GetXaxis()->GetBinLowEdge(ix+1) + angular2DResponse_hits[i][p][ibinx]->GetXaxis()->GetBinWidth(ix+1);
	      double yposition = angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinCenter(iy+1);
	      if(iy == 0) yposition = angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinLowEdge(iy+1);
	      else if(iy == (angular2DResponse_hits[i][p][ibinx]->GetNbinsY()-1)) yposition = angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinLowEdge(iy+1) + angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinWidth(iy+1);
	      gAngular2DResponse[i][p][ibinx]->SetPoint(ig,xposition,yposition,content);
	      ig++;
	      //gAngular2DResponse[i][p][ibinx]->SetPoint(ig,angular2DResponse_hits[i][p][ibinx]->GetXaxis()->GetBinLowEdge(ix+1),angular2DResponse_hits[i][p][ibinx]->GetYaxis()->GetBinLowEdge(iy+1),content);
	    }
	  }

	  //1. Normalize
	  angular2DResponse_hits[i][p][ibinx]->Scale(1./angular2DResponse_hits[i][p][ibinx]->Integral());

	  //gAngular2DResponse[i][p][ibinx]->Draw("colz");
	  angular2DResponse_hits[i][p][ibinx]->Draw("colz");
	  cAngularResponsePMT_hits[i][p][ibinx]->SaveAs(Form("plots/cAngularResponsePMT2D_pmtType%d_group%d_distance%d.eps",i,p,ibinx));
	  angular2DResponse_hits[i][p][ibinx]->Write(Form("hPMTDirectionality_2D_%d_%d_%d",0,i,p));
	  gAngular2DResponse[i][p][ibinx]->Write(Form("gPMTDirectionality_2D_%d_%d_%d",0,i,p));
	  
	  //angular2DResponse_hits[i][p][ibinx]->Write(Form("angularPDF_pmtType%d_group%d",i,p));
	  //delete min;
	}



	


	
      }
    }
 
    if(f==0) cChargePerPMT = new TCanvas("cChargePerPMT","cChargePerPMT");
    for(int i=0;i<nPMTtypes;i++){
      ChargePerPMT[f][i]->Scale(1./ChargePerPMT[f][i]->Integral());
      if(i==0){
	ChargePerPMT[f][i]->Draw("hist");
	ChargePerPMT[f][i]->GetXaxis()->SetRangeUser(0,8);
	ChargePerPMT[f][i]->GetYaxis()->SetRangeUser(0,1.);
      }
      else  ChargePerPMT[f][i]->Draw("histsame");
      l->Draw("same");
    }
    cChargePerPMT->SaveAs(Form("plots/%s.eps",cChargePerPMT->GetName()));

    
    if(f==0){
      cChargeProfile = new TCanvas("cChargeProfile","cChargeProfile");
      cChargeProfile->Divide(1,2);
    }
    for(int i=0;i<nPMTtypes;i++){
      ChargeProfile[f][i]->Scale(1/nEvents);
      ChargeProfile[f][i]->Rebin(4);

      ratioChargeProfile[f][i] = (TH1D*) ChargeProfile[f][i]->Clone(Form("ratioChargeProfile_%d_%d",f,i));
      ratioChargeProfile[f][i]->Divide(ChargeProfile[f][0]);

      if(i==0){
	cChargeProfile->cd(1);
	ChargeProfile[f][i]->Draw("hist");
	//ChargeProfile[f][i]->GetXaxis()->SetRangeUser(0,8);
	//ChargeProfile[f][i]->GetYaxis()->SetRangeUser(0,1.);
      }
      else{
	cChargeProfile->cd(1);
	ChargeProfile[f][i]->Draw("histsame");
	l->Draw("same");
	if(i==1){
	  cChargeProfile->cd(2);
	  ratioChargeProfile[f][i]->Draw("hist");
	}
	else{
	  cChargeProfile->cd(2);
	  ratioChargeProfile[f][i]->Draw("histsame");
	  l->Draw("same");
	}
      }
    }
    cChargeProfile->SaveAs(Form("plots/%s.eps",cChargeProfile->GetName()));


    if(isNew2){
      if(f==0){
	cChargeProfileXdWall = new TCanvas("cChargeProfileXdWall","cChargeProfileXdWall");
	cChargeProfileXdWall->Divide(1,2);
      }
      for(int i=0;i<nPMTtypes;i++){
	ChargeProfileXdWall[f][i]->Scale(1/nEvents);
	ChargeProfileXdWall[f][i]->Rebin2D(100,4);
	
	cChargeProfileXdWall->cd(i+1);
	ChargeProfileXdWall[f][i]->Draw("colz");
      }
      cChargeProfileXdWall->SaveAs(Form("plots/%s.eps",cChargeProfileXdWall->GetName()));
    }
    
    if(f==0) cTotalCharge = new TCanvas("cTotalCharge","cTotalCharge");
    for(int i=0;i<nPMTtypes;i++){
      TotalCharge[f][i]->Scale(1./TotalCharge[f][i]->Integral());
      if(i==0){
	TotalCharge[f][i]->Draw("hist");
	TotalCharge[f][i]->GetXaxis()->SetRangeUser(0,5e3);
	TotalCharge[f][i]->GetYaxis()->SetRangeUser(0,0.6);
      }
      else  TotalCharge[f][i]->Draw("histsame");
      l->Draw("same");
    }
    cTotalCharge->SaveAs(Form("plots/%s.eps",cTotalCharge->GetName()));

    
    if(f==0) cTotalHit = new TCanvas("cTotalHit","cTotalHit");
    for(int i=0;i<nPMTtypes;i++){
      TotalHit[f][i]->Scale(1./TotalHit[f][i]->Integral());
      if(i==0){
	TotalHit[f][i]->Draw("hist");
	TotalHit[f][i]->GetXaxis()->SetRangeUser(0,5e3);
	TotalHit[f][i]->GetYaxis()->SetRangeUser(0,0.6);
      }
      else  TotalHit[f][i]->Draw("histsame");
      l->Draw("same");
    }
    cTotalHit->SaveAs(Form("plots/%s.eps",cTotalHit->GetName()));


    if(f==0) cTimeProfile = new TCanvas("cTimeProfile","cTimeProfile");
    for(int i=0;i<nPMTtypes;i++){
      TimeProfile[f][i]->Scale(1/TimeProfile[f][i]->Integral());
      //TimeProfile[f][i]->Scale(1/nEvents);
      TimeProfile[f][i]->Rebin(4);
      if(i==0){
	TimeProfile[f][i]->Draw("hist");
	TimeProfile[f][i]->GetXaxis()->SetRangeUser(0,600);
	//TimeProfile[f][i]->GetYaxis()->SetRangeUser(0,1.);
      }
      else  TimeProfile[f][i]->Draw("histsame");
      l->Draw("same");
    }
    cTimeProfile->SaveAs(Form("plots/%s.eps",cTimeProfile->GetName()));



    //Very important here: let profile and DR have same transformation (same binning, same rescaling etc) as they should be kept in the same proportions as signal
    if(f==0) cTimeTOFProfile = new TCanvas("cTimeTOFProfile","cTimeTOFProfile");
    for(int i=0;i<nPMTtypes;i++){
      //We will zoom in the time tof profile region from -100 to -20ns.
      
      TimeTOFProfile[f][i]->Rebin(1);
      //cout<<"Integral signal from -100 to 500ns = "<<TimeTOFProfile[f][i]->Integral(TimeTOFProfile[f][i]->FindBin(-100),TimeTOFProfile[f][i]->FindBin(500))<<", DR="<<TimeTOFDR[f][i]->Integral(TimeTOFDR[f][i]->FindBin(-100),TimeTOFDR[f][i]->FindBin(500))<<endl;
    
      double scaleValue=events[i]/2;//TimeTOFProfile[f][i]->GetBinContent(TimeTOFProfile[f][i]->GetMaximumBin());
      TimeTOFProfile[f][i]->Scale(1/scaleValue);
      if(i==0){
	TimeTOFProfile[f][i]->Draw("hist");
	TimeTOFProfile[f][i]->GetXaxis()->SetRangeUser(-10,30);
	//TimeTOFProfile[f][i]->GetYaxis()->SetRangeUser(0,1.);
      }
      else TimeTOFProfile[f][i]->Draw("histsame");
    }   
    ///////////////////////////////////////
  
      //Very important here: let profile and DR have same transformation (same binning, same rescaling etc) as they should be kept in the same proportions as signal
      if(f==0) cHitTimeTOFProfile = new TCanvas("cHitTimeTOFProfile","cHitTimeTOFProfile");
      for(int i=0;i<nPMTtypes;i++){
	//We will zoom in the time tof profile region from -100 to -20ns.
	double startDR=-100;
	double endDR=-20;
	double timeWindowDR=endDR-startDR;
	DRTotalPerNS[i]=HitTimeTOFProfile[f][i]->Integral(HitTimeTOFProfile[f][i]->FindBin(startDR),HitTimeTOFProfile[f][i]->FindBin(endDR));
	DRTotalPerNS[i]/=timeWindowDR;
	
	for(int ibinx=1;ibinx<= HitTimeTOFProfile[f][i]->GetNbinsX();ibinx++){
	  double timeWindow=HitTimeTOFProfile[f][i]->GetBinWidth(ibinx);
	  HitTimeTOFDR[f][i]->SetBinContent(ibinx,DRTotalPerNS[i]*timeWindow);
	  //cout<<HitTimeTOFDR[f][i]->GetBinCenter(ibinx)<<", time window="<<timeWindow<<", signal="<<HitTimeTOFProfile[f][i]->GetBinContent(ibinx)<<", DR="<<HitTimeTOFDR[f][i]->GetBinContent(ibinx)<<", bin width="<<timeWindow<<", DR="<<DRTotalPerNS[i]<<", nevents="<<events[i]/2<<endl;
	}
	HitTimeTOFProfile[f][i]->Rebin(1);
	HitTimeTOFDR[f][i]->Rebin(1);
	//cout<<"Integral signal from -100 to 500ns = "<<HitTimeTOFProfile[f][i]->Integral(HitTimeTOFProfile[f][i]->FindBin(-100),HitTimeTOFProfile[f][i]->FindBin(500))<<", DR="<<HitTimeTOFDR[f][i]->Integral(HitTimeTOFDR[f][i]->FindBin(-100),HitTimeTOFDR[f][i]->FindBin(500))<<endl;
	
	double scaleValue=events[i]/2;//HitTimeTOFProfile[f][i]->GetBinContent(HitTimeTOFProfile[f][i]->GetMaximumBin());
	HitTimeTOFProfile[f][i]->Scale(1/scaleValue);
	HitTimeTOFDR[f][i]->Scale(1/scaleValue);
	//HitTimeTOFProfile[f][i]->Scale(1/nEvents);
	if(i==0){
	  HitTimeTOFProfile[f][i]->Draw("hist");
	  HitTimeTOFProfile[f][i]->GetXaxis()->SetRangeUser(-10,30);
	//HitTimeTOFProfile[f][i]->GetYaxis()->SetRangeUser(0,1.);
	}
      else HitTimeTOFProfile[f][i]->Draw("histsame");
	HitTimeTOFDR[f][i]->SetLineColor(kMagenta);
	HitTimeTOFDR[f][i]->Draw("histsame");
	cout<<"Drawn"<<endl;
	cout<<"Integral signal from -100 to 500ns = "<<HitTimeTOFProfile[f][i]->Integral(HitTimeTOFProfile[f][i]->FindBin(-100),HitTimeTOFProfile[f][i]->FindBin(500))<<", DR="<<HitTimeTOFDR[f][i]->Integral(HitTimeTOFDR[f][i]->FindBin(-100),HitTimeTOFDR[f][i]->FindBin(500))<<endl;
   
      ///////////////////////////////////////
      
      gausExpoConv[f][i] = new TF1(Form("gausExpoConv%d_%d",f,i),".5*[0]*TMath::Exp(-[3]*((x-[1])-[2]*[2]*[3]/2))*(1 + TMath::Erf( ((x-[1])-[2]*[2]*[3])/(TMath::Sqrt(2)*[2])))",-2,limitFitGausExpo[i]);
      gausExpoConv[f][i]->SetLineColor(1);
      gausExpoConv[f][i]->SetLineWidth(2);
      gausExpoConv[f][i]->SetParameter(0,10.);
      gausExpoConv[f][i]->SetParameter(1,0.);
      gausExpoConv[f][i]->SetParameter(2,2);
      gausExpoConv[f][i]->SetParameter(3,1.);
      HitTimeTOFProfile[f][i]->Fit(Form("gausExpoConv%d_%d",f,i),"R");
      gausExpoConv[f][i]->Draw("same");
      l->Draw("same");
      
      // int nBins=
      //double xPos
      //for(double 
      //TGraph * graphsplineExpoConv[f][i] = new TGraph(Form("graphExpoConv%d_%d",f,i),nBins,xPos,yPos);

      splineExpoConv[f][i] = new TSpline3(Form("splineExpoConv%d_%d",f,i),-limitFitGausExpo[i],limitFitGausExpo[i],gausExpoConv[f][i],300);
      splineExpoConv[f][i]->SetLineColor(kCyan);
      //splineExpoConv[f][i]->Draw("lcsame");
      fOut->cd();
      splineExpoConv[f][i]->Write(splineExpoConv[f][i]->GetTitle());

      cout<<"Start to merge bins in higher size bins to optimize PDF"<<endl;
      vector <double> xPos;xPos.clear();
      vector <double> yPos;yPos.clear();
      vector <double> drPos;drPos.clear();
      int iBinActive=0;
      int nRebins=10;//Number of small bin gathered in a big one
      double xAverage=0;
      double yAverage=0;
      double drAverage=0;
      int binLimit=HitTimeTOFProfile[f][i]->FindBin(-limitFitGausExpo[i]);//HitTimeTOFProfile[f][i]->FindBin(limitFitGausExpo[i]-5);//Lower limit of the active big bin
      const int nLimits = 8;
      double minValue = -limitFitGausExpo[i];
      double maxValue = HitTimeTOFProfile[f][i]->GetXaxis()->GetXmax();
      double lowLimits[nLimits] = {minValue, -100., -30., -10., 8., 30., 100.,maxValue};
      double binLowLimits[nLimits];
      for(int il = 0;il<nLimits;il++){
	binLowLimits[il] = HitTimeTOFProfile[f][i]->FindBin(lowLimits[il]);
      }
      int rebinLowLimits[nLimits] = {100, 10, 5, 2, 10, 30, 100, 100};
      int currentLimit = 0;
      int nBinsAverage = 0;
      int nBins=0;//150;//Number of big bins
      //int iBinActive = 0;
      //In between the different limits, the rebinning factor will be different.
      //We should first ensure in which region we are. Based on this, we use a different rebinning value.
      //Then, we will loop over the bins which are in the correct region. We will have an internal counter. As soon as this internal counter reach the rebinning value, we stop and average.
      //Another way of reaching the factor could be to reach another regiob. We should take this into account.
	
      //for(int ibinx=HitTimeTOFProfile[f][i]->FindBin(limitFitGausExpo[i]-5);ibinx<=HitTimeTOFProfile[f][i]->GetNbinsX();ibinx++){
      //Loop over all bins of the histogram
      for(int ibinx=HitTimeTOFProfile[f][i]->FindBin(-limitFitGausExpo[i]);ibinx<=HitTimeTOFProfile[f][i]->GetNbinsX();ibinx++){
	cout<<"Bin = "<<ibinx<<" i.e. value of low edge = "<< HitTimeTOFProfile[f][i]->GetBinLowEdge(ibinx)<<", current low edge limit bin = "<<binLowLimits[currentLimit]<<", nBinsAverage = "<<nBinsAverage<<", rebinning factor = "<<rebinLowLimits[currentLimit]<<endl;
	if((ibinx >= binLowLimits[currentLimit] && ibinx < binLowLimits[currentLimit+1]) && nBinsAverage < rebinLowLimits[currentLimit]){
	  //Sum over bin content to average them until we reach the last bin 
	  xAverage+=HitTimeTOFProfile[f][i]->GetBinCenter(ibinx);
	  yAverage+=HitTimeTOFProfile[f][i]->GetBinContent(ibinx);
	  drAverage+=HitTimeTOFDR[f][i]->GetBinContent(ibinx);
	  nBinsAverage++;
	}
	if( (ibinx+1 >= binLowLimits[currentLimit+1]) || (nBinsAverage >= rebinLowLimits[currentLimit]) ){
	  //If we reached the number of binning to rebin: store information.
	  //Same if we reached the limit of the region to rebin of a given factor.
	  //What if both happens at the same time, or worse: we reach limit of bins on bin n, and at n+1, we overcome the limit. In that case, we do not have anyting to fill our average? So, we understand that if we are in the last bin below the limit, we should store and then pass to the next step of limit. So, the check of the limit should always be on the next bin. 
	  xPos.push_back(xAverage/nBinsAverage);
	  drPos.push_back(drAverage/nBinsAverage);
	  yPos.push_back((yAverage/nBinsAverage));
	  cout<<"Position = "<<xPos[nBins]<<", value="<<yPos[nBins]<<", DR = "<<drPos[nBins]<<endl;
	  nBins++;
	  nBinsAverage = 0;
	  xAverage=0;
	  yAverage=0;
	  drAverage=0;
	  if(ibinx+1 >= binLowLimits[currentLimit+1]){
	    currentLimit++;
	  }
	}
      }
      double nBins_graph = nBins;
      double * xPos_graph = new double[nBins];
      double * yPos_graph = new double[nBins];
      double * drPos_graph = new double[nBins];
      for(int ibin = 0; ibin<nBins;ibin++){
	xPos_graph[ibin] = xPos.at(ibin);
	yPos_graph[ibin] = yPos.at(ibin);
	drPos_graph[ibin] = drPos.at(ibin);
      }
      graphExpoQueue[f][i] = new TGraph(nBins_graph,xPos_graph,yPos_graph);
      
      //HitTimeTOFProfile[f][i]->Rebin(4);
      //ExpoQueue[f][i] = new TF1(Form("ExpoQueue%d_%d",f,i),"[0]*TMath::Exp(-[1]*x)",limitFitGausExpo[i],500);
      //ExpoQueue[f][i]->SetLineColor(kGray);
      //ExpoQueue[f][i]->SetLineWidth(2);
      //ExpoQueue[f][i]->SetParameter(0, HitTimeTOFProfile[f][i]->GetBinContent(HitTimeTOFProfile[f][i]->FindBin(limitFitGausExpo[i])));
      //ExpoQueue[f][i]->SetParameter(1,1.);
      //HitTimeTOFProfile[f][i]->Fit(Form("ExpoQueue%d_%d",f,i),"R");
      //ExpoQueue[f][i]->Draw("same");

      //splineExpoQueue[f][i] = new TSpline3(Form("splineExpoQueue%d_%d",f,i),limitFitGausExpo[i],500.,ExpoQueue[f][i],300);
      splineExpoQueue[f][i] = new TSpline3(Form("splineExpoQueue%d_%d",f,i),graphExpoQueue[f][i]);
      splineExpoQueue[f][i]->SetLineColor(kCyan);
      graphExpoQueue[f][i]->SetMarkerColor(kCyan);
      graphExpoQueue[f][i]->SetMarkerStyle(20);
      graphExpoQueue[f][i]->Draw("Psame");
      splineExpoQueue[f][i]->Draw("lcsame");
      fOut->cd();
      splineExpoQueue[f][i]->Write(splineExpoQueue[f][i]->GetTitle());


      graphDR[f][i] = new TGraph(nBins_graph,xPos_graph,drPos_graph);
      //splineDR[f][i] = new TSpline3(Form("splineDR%d_%d",f,i),limitFitGausExpo[i],500.,DR[f][i],300);
      splineDR[f][i] = new TSpline3(Form("splineDR%d_%d",f,i),graphDR[f][i]);
      splineDR[f][i]->SetLineColor(kCyan);
      splineDR[f][i]->Draw("lcsame");
      fOut->cd();
      splineDR[f][i]->Write(splineDR[f][i]->GetTitle());
      
      cout<<"Integral signal from -100 to 500ns = "<<HitTimeTOFProfile[f][i]->Integral(HitTimeTOFProfile[f][i]->FindBin(-100),HitTimeTOFProfile[f][i]->FindBin(500))<<", DR="<<HitTimeTOFDR[f][i]->Integral(HitTimeTOFDR[f][i]->FindBin(-100),HitTimeTOFDR[f][i]->FindBin(500))<<endl;
      double stepSize=0.5;
      double integralSignal=0;
      double integralDR=0;
      for(double j=-100;j<500;j+=stepSize){
	integralSignal+=stepSize*splineExpoQueue[f][i]->Eval(j+stepSize/2);
	integralDR+=stepSize*splineDR[f][i]->Eval(j+stepSize/2);
      }
      cout<<"Integral spline="<<integralSignal<<", DR="<<integralDR<<endl;
      }
    cHitTimeTOFProfile->SaveAs(Form("plots/%s.eps",cHitTimeTOFProfile->GetName()));
    
    if(f==0){
      cTimeTOFProfileXTOF = new TCanvas("cTimeTOFProfileXTOF","cTimeTOFProfileXTOF");
      cTimeTOFProfileXTOF->Divide(1,2);
    
      for(int i=0;i<nPMTtypes;i++){
	TimeTOFProfileXTOF[f][i]->Scale(1/TimeTOFProfileXTOF[f][i]->Integral());
	//TimeTOFProfileXTOF[f][i]->Scale(1/nEvents);
	TimeTOFProfileXTOF[f][i]->Rebin2D(2,4);
	//TimeTOFProfileXTOF[f][i]->GetYaxis()->Rebin(3);
	cTimeTOFProfileXTOF->cd(i+1);
	TimeTOFProfileXTOF[f][i]->Draw("colz");
	TimeTOFProfileXTOF[f][i]->GetXaxis()->SetRangeUser(-50,100);
	TimeTOFProfileXTOF[f][i]->GetYaxis()->SetRangeUser(0,300);
      }
      cTimeTOFProfileXTOF->SaveAs(Form("plots/%s.eps",cTimeTOFProfileXTOF->GetName()));
    }

        
    if(f==0){
      cTimeTOFProfileXTOF_1D = new TCanvas("cTimeTOFProfileXTOF_1D","cTimeTOFProfileXTOF_1D");
      //cTimeTOFProfileXTOF_1D->Divide( floor(TMath::Sqrt(nBinsTOF)) , floor(TMath::Sqrt(nBinsTOF)) );
      cTimeTOFProfileXTOF_1D->Divide( floor(TMath::Sqrt(nBinsTOF)) , floor(TMath::Sqrt(nBinsTOF)) );
    
      for(int i=0;i<nPMTtypes;i++){
	
	for(int ibiny=1;ibiny<=nBinsTOF;ibiny++){

	  TimeTOFProfileXTOF_1D[f][i][ibiny-1] = (TH1D*) TimeTOFProfileXTOF[f][i]->ProjectionX(Form("TimeTOFProfileXTOF%d_%d_%d",f,i,ibiny),ibiny,ibiny);
	  //if(ibiny<26) continue;
	  //else
	  cTimeTOFProfileXTOF_1D->cd(ibiny);
	  TimeTOFProfileXTOF_1D[f][i][ibiny-1]->SetLineColor(colorScale[i]);

	  if(i==0){
	    TimeTOFProfileXTOF_1D[f][i][ibiny-1]->Draw("hist");
	    TimeTOFProfileXTOF_1D[f][i][ibiny-1]->GetXaxis()->SetRangeUser(-5,10);
	    //TimeTOFProfile[f][i]->GetYaxis()->SetRangeUser(0,1.);
	  }
	  else  TimeTOFProfileXTOF_1D[f][i][ibiny-1]->Draw("histsame");
	  l->Draw("same");

	  gausExpoConv_slice[f][i][ibiny-1] = new TF1(Form("gausExpoConv_slice%d_%d_%d",f,i,ibiny-1),".5*[0]*TMath::Exp(-[3]*((x-[1])-[2]*[2]*[3]/2))*(1 + TMath::Erf( ((x-[1])-[2]*[2]*[3])/(TMath::Sqrt(2)*[2])))",-3,10);
	  gausExpoConv_slice[f][i][ibiny-1]->SetLineColor(1);
	  gausExpoConv_slice[f][i][ibiny-1]->SetLineWidth(2);
	  gausExpoConv_slice[f][i][ibiny-1]->SetParameter(0,10.);
	  gausExpoConv_slice[f][i][ibiny-1]->SetParameter(1,0.);
	  gausExpoConv_slice[f][i][ibiny-1]->SetParameter(2,1.);
	  gausExpoConv_slice[f][i][ibiny-1]->SetParameter(3,2.);

	  TimeTOFProfileXTOF_1D[f][i][ibiny-1]->Fit(Form("gausExpoConv_slice%d_%d_%d",f,i,ibiny-1),"R");
	  gausExpoConv_slice[f][i][ibiny-1]->Draw("same");
	}
      }
    }
  }
    
  fOut->Close();
}

