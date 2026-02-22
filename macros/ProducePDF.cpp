#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TChain.h>
#include <TBranch.h>
#include <TF1.h>
#include <TSpline.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMath.h>


using namespace std;




struct arguments
{
  char * inputFile  = NULL;
  char * outputFile = NULL;
  bool verbose      = false;
};

arguments FetchInput(int argc, char* argv[])
{
    int c = -1;
    arguments arglist;

    //Input in c the argument (-f etc...) and in optarg the next argument
    //When the above test becomes -1, it means it fails to find a new argument
    std::cout << std::endl;
    while( (c = getopt(argc, argv, "i:o:v")) != -1 )
    {
        switch(c)
        {
        //Input file name
        case 'i':
            arglist.inputFile = optarg;
            std::cout << "Input WCSim file: " << arglist.inputFile << std::endl;
            break; 

        //Output file name
        case 'o':
            arglist.outputFile = optarg;
            if(arglist.outputFile == NULL){sprintf(arglist.outputFile,"PDF.root");}
            std::cout << "Output root file: " << arglist.outputFile << std::endl;
            break;

        //Warning
        case 'v':
            arglist.verbose = true;
            std::cout << "VERBOSE option on" << std::endl;
            break;
        }
    }
    std::cout << std::endl;

    return arglist;
}



/*****************************************************************************************************/
/* MAKE PDFs */
/*****************************************************************************************************/

int main(int argc, char **argv)
{
  const int nFiles=1;
  TFile * _f;

  //Get arguments
  arguments arglist = FetchInput(argc, argv);
  char * ifname = arglist.inputFile;
  char * ofname = arglist.outputFile;
  bool verbose  = arglist.verbose;

  if (ifname==NULL)
  {
    cout << "Error, no input file" << endl;
    return -1;
  }
  else
  {
    _f = new TFile(ifname,"read");
  }

  if (!_f->IsOpen())
  {
    cout << "Error, could not open input file: " << ifname << endl;
    return -1;
  }

  TFile * fOut = new TFile(ofname,"recreate");

  int opt;
  const int nPMTtypes = 2;
  double DRTotalPerNS[nPMTtypes];
  double limitFitGausExpo[nPMTtypes];
  for(int i=0;i<nPMTtypes;i++)
  {
    if(i==0) limitFitGausExpo[i] = 800;
    else limitFitGausExpo[i]=800;
  }

  // Define TH1D, TGraph and TSpline for TOF residual PDF
  TH1D * TotalCharge[nPMTtypes];
  TH1D * TimeTOFProfile[nPMTtypes];
  TH1D * HitTimeTOFProfile[nPMTtypes];
  TH1D * HitTimeTOFDR[nPMTtypes];
  
  TGraph * graphExpoConv[nPMTtypes];
  TSpline3 * splineExpoConv[nPMTtypes];

  TGraph * graphExpoQueue[nPMTtypes];
  TSpline3 * splineExpoQueue[nPMTtypes];

  TGraph * graphDR[nPMTtypes];
  TSpline3 * splineDR[nPMTtypes];

  // Define TH1D, TGraph and TSpline for Cherenkov angle PDF
  TH1D * ChargeProfile[nPMTtypes];
  TGraph * graphChargeProfile[nPMTtypes];
  TSpline3 * splineChargeProfile[nPMTtypes];


  for(int i=0 ; i<nPMTtypes ; i++)
  {
    ChargeProfile[i]     = (TH1D*) _f->Get(Form("ChargeProfile_pmtType%d", i));
    TotalCharge[i]       = (TH1D*) _f->Get(Form("TotalCharge_pmtType%d", i));
    HitTimeTOFProfile[i] = (TH1D*) _f->Get(Form("HitTimeTOFProfile_pmtType%d", i));
    HitTimeTOFDR[i]      = (TH1D*) HitTimeTOFProfile[i]->Clone(Form("HitTimeTOFDR_pmtType%d", i));
    HitTimeTOFDR[i]->Reset();   // We rest because we only care about the TH1D structure and not really about the content since we recompute the content below

    //Very important here: let profile and DR have same transformation (same binning, same rescaling etc) as they should be kept in the same proportions as signal
    //We will zoom in the time tof profile region from -100 to -20ns.
    double startDR = -100;
    double endDR   = -20;
    double timeWindowDR = endDR - startDR;
    DRTotalPerNS[i]  = HitTimeTOFProfile[i]->Integral(HitTimeTOFProfile[i]->FindBin(startDR), HitTimeTOFProfile[i]->FindBin(endDR));
    DRTotalPerNS[i] /= timeWindowDR;
    
    // Set dark rate to the average dark rate in [-100, -20] ns
    for(int ibinx=1;ibinx<= HitTimeTOFProfile[i]->GetNbinsX();ibinx++)
    {
      double timeWindow=HitTimeTOFProfile[i]->GetBinWidth(ibinx);
      HitTimeTOFDR[i]->SetBinContent(ibinx, DRTotalPerNS[i]*timeWindow);
    }
    
    // Scaling with arbitrary factor to get an approximate PDF for HitTimeTOFProfile
    // HitTimeTOFDR is scaled by the same factor to keep both PDF consistent with each other
    double scaleValue = TotalCharge[i]->Integral()/2;
    HitTimeTOFProfile[i]->Scale(1/scaleValue);
    HitTimeTOFDR[i]->Scale(1/scaleValue);


    //////////////////////////////////////////////
    // MAKE SPLINE OVER ALL TOF RESIDUALS RANGE //
    //////////////////////////////////////////////

    int nBinsConv = HitTimeTOFProfile[i]->GetNbinsX();
    double * xConv_graph  = new double[nBinsConv];
    double * yConv_graph  = new double[nBinsConv];
    for(int ibin = 0; ibin<nBinsConv ; ibin++)
    {
      xConv_graph[ibin] = HitTimeTOFProfile[i]->GetBinCenter(ibin);
      yConv_graph[ibin] = HitTimeTOFProfile[i]->GetBinContent(ibin);
    }

    graphExpoConv[i]  = new TGraph(nBinsConv, xConv_graph, yConv_graph);
    splineExpoConv[i] = new TSpline3(Form("splineExpoConv_%d", i), graphExpoConv[i]);
    fOut->cd();
    splineExpoConv[i]->Write(splineExpoConv[i]->GetTitle());


    ////////////////
    // MERGE BINS //
    ////////////////

    cout << "Merge bins to optimize PDF" << endl;
    vector <double> xPos; xPos.clear();
    vector <double> yPos; yPos.clear();
    vector <double> drPos;drPos.clear();
    int iBinActive=0;
    int nRebins=10;   //Number of small bin gathered in a big one
    double xAverage=0;
    double yAverage=0;
    double drAverage=0;
    int binLimit=HitTimeTOFProfile[i]->FindBin(-limitFitGausExpo[i]);   //Lower limit of the active big bin
    const int nLimits = 8;
    double minValue = -limitFitGausExpo[i];
    double maxValue = HitTimeTOFProfile[i]->GetXaxis()->GetXmax();
    double lowLimits[nLimits] = {minValue, -100., -30., -10., 8., 30., 100.,maxValue};
    double binLowLimits[nLimits];
    for(int il=0 ; il<nLimits ; il++){binLowLimits[il] = HitTimeTOFProfile[i]->FindBin(lowLimits[il]);}
    int rebinLowLimits[nLimits] = {100, 10, 5, 2, 10, 30, 100, 100};
    int currentLimit = 0;
    int nBinsAverage = 0;
    int nBins=0;   //Number of bigger bins after merging

    // In between the different limits, the rebinning factor will be different.
    // We should first ensure in which region we are. Based on this, we use a different rebinning value.
    // Then, we will loop over the bins which are in the correct region. We will have an internal counter.
    // As soon as this internal counter reach the rebinning value, we stop and average.
    // Another way of reaching the factor could be to reach another region. We should take this into account.

    for(int ibinx=HitTimeTOFProfile[i]->FindBin(-limitFitGausExpo[i]) ; ibinx<=HitTimeTOFProfile[i]->GetNbinsX() ; ibinx++)
    {
      cout << "Bin = " << ibinx
           << " i.e. value of low edge = "<< HitTimeTOFProfile[i]->GetBinLowEdge(ibinx)
           << ", current low edge limit bin = " << binLowLimits[currentLimit]
           << ", nBinsAverage = " << nBinsAverage
           << ", rebinning factor = " << rebinLowLimits[currentLimit] << endl;

      // Sum over bin content to average them until we reach the last bin 
      if((ibinx >= binLowLimits[currentLimit] && ibinx < binLowLimits[currentLimit+1]) && nBinsAverage < rebinLowLimits[currentLimit])
      {
        xAverage  += HitTimeTOFProfile[i]->GetBinCenter(ibinx);
        yAverage  += HitTimeTOFProfile[i]->GetBinContent(ibinx);
        drAverage += HitTimeTOFDR[i]->GetBinContent(ibinx);
        nBinsAverage++;
      }

      // If we have reached the number of bins to rebin: store information.
      // Same if we have reached the limit of the region to rebin of a given factor.
      // What if both happens at the same time, or worse: we reach limit of bins on bin n, and at n+1, we overcome the limit. In that case, we do not have anyting to fill our average?
      // So, we understand that if we are in the last bin below the limit, we should store and then pass to the next step of limit. So, the check of the limit should always be on the next bin.

      if( (ibinx+1 >= binLowLimits[currentLimit+1]) || (nBinsAverage >= rebinLowLimits[currentLimit]) )
      {
        xPos.push_back( xAverage/nBinsAverage);
        yPos.push_back( yAverage/nBinsAverage);
        drPos.push_back(drAverage/nBinsAverage);

        cout << "Position = " << xPos[nBins] << ", value=" << yPos[nBins] << ", DR = " << drPos[nBins] << endl;

        nBins++;
        nBinsAverage = 0;
        xAverage     = 0;
        yAverage     = 0;
        drAverage    = 0;
        if(ibinx+1 >= binLowLimits[currentLimit+1]){currentLimit++;}
      }
    }

    double * xPos_graph  = new double[nBins];
    double * yPos_graph  = new double[nBins];
    double * drPos_graph = new double[nBins];
    for(int ibin = 0; ibin<nBins ; ibin++)
    {
      xPos_graph[ibin]  = xPos.at(ibin);
      yPos_graph[ibin]  = yPos.at(ibin);
      drPos_graph[ibin] = drPos.at(ibin);
    }


    //////////////////////////////////
    // MAKE SPLINE FOR QUEUE SIGNAL //
    //////////////////////////////////

    graphExpoQueue[i]  = new TGraph(nBins, xPos_graph, yPos_graph);
    splineExpoQueue[i] = new TSpline3(Form("splineExpoQueue_%d", i), graphExpoQueue[i]);
    fOut->cd();
    splineExpoQueue[i]->Write(splineExpoQueue[i]->GetTitle());


    ////////////////////////////////
    // MAKE SPLINE FOR DARK NOISE //
    ////////////////////////////////

    graphDR[i] = new TGraph(nBins, xPos_graph, drPos_graph);
    splineDR[i] = new TSpline3(Form("splineDR_%d", i), graphDR[i]);
    fOut->cd();
    splineDR[i]->Write(splineDR[i]->GetTitle());


    /////////////////////////////////////
    // MAKE SPLINE FOR CHERENKOV ANGLE //
    /////////////////////////////////////

    int nBinsAngle = ChargeProfile[i]->GetNbinsX();
    double * xAngle_graph  = new double[nBinsAngle];
    double * yAngle_graph  = new double[nBinsAngle];
    for(int ibin = 0; ibin<nBinsAngle ; ibin++)
    {
      xAngle_graph[ibin] = ChargeProfile[i]->GetBinCenter(ibin);
      yAngle_graph[ibin] = ChargeProfile[i]->GetBinContent(ibin);
    }

    graphChargeProfile[i]  = new TGraph(nBinsAngle, xAngle_graph, yAngle_graph);
    splineChargeProfile[i] = new TSpline3(Form("splineChargeProfile_%d", i), graphChargeProfile[i]);
    fOut->cd();
    splineChargeProfile[i]->Write(splineChargeProfile[i]->GetTitle());
  }
  
  fOut->Close();
}

