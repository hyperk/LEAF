// This is a sample script to be executed in root
// Usage: root generateShell_analyzeWCSim.c

{
  char * ShellName = new char[256];
  ofstream ShellFile(ShellName);
  char * infilename = new char[512];char * infilename2 = new char[512];
  char * outfilename = new char[512];char * outfilename2 = new char[512];
  char * LEAFDIR = getenv("LEAFDIR");

  int nevents=100000;
  int neventsPerFile=1000;

  for(int nrj=10;nrj<11;nrj++){
    sprintf(infilename,"/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_fulltank_0hitstrigger_100000_nodr.root",nrj); // WCSim file to be analyzed

    for(int startBin=0;startBin<nevents;startBin+=neventsPerFile){
      sprintf(ShellName,"%s/jobs/analyzeLE_e%d_start%d.sh",LEAFDIR,nrj,startBin); // Shell file to be submitted in queue system
      sprintf(outfilename,"/disk01/usr5/bquilain/Analyze_hkhybridmpmt10pc14374100Hz_4.2kHzbl_%dMeV_tmp%d.root",nrj,startBin);

      ofstream ShellFile(ShellName);
      if(ShellFile){
        ShellFile<<"#!/bin/bash"<<endl;
        ShellFile<<Form("rm -f %s",outfilename)<<endl;
        ShellFile<<Form("%s/macros/AnalyzeWSHierarchy -h -f %s -o %s -s %d -e %d",LEAFDIR,infilename,outfilename,startBin,startBin+neventsPerFile)<<endl;
      }
    }
  }
}
