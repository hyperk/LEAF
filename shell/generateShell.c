// This is a sample script to be executed in root
// Usage: root generateShell.c

{
  char * ShellName = new char[256];
  ofstream ShellFile(ShellName);
  char * infilename = new char[512];char * infilename2 = new char[512];
  char * outfilename = new char[512];char * outfilename2 = new char[512];
  char * LEAFDIR = getenv("LEAFDIR");
  double dr = 0;
  bool mPMT = true;
  double dr_mPMT = 0;

  int nevents=50000;
  int neventsPerFile=500;
  for(int nrj=3;nrj<16;nrj++){ // nrj -> energy
    sprintf(infilename,"/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_fulltank_0hitstrigger_100000.root",nrj); // WCSim file to be analyzed

    for(int startBin=0;startBin<nevents;startBin+=neventsPerFile){
      sprintf(ShellName,"%s/jobs/fitLE_e%d_start%d.sh",LEAFDIR,nrj,startBin); // shell file to be submitted to queue system
      sprintf(outfilename,"%s/jobs/root/hkhybridmpmt10pc14374100Hz_4.2kHzbl_%dMeV_tmp%d.root",LEAFDIR, nrj,startBin); // outpu t file name of process 

      ofstream ShellFile(ShellName);
      if(ShellFile){
        ShellFile<<"#!/bin/bash"<<endl;
        ShellFile<<Form("rm -f %s",outfilename)<<endl;
        ShellFile<<Form("source %s/RunAtStart_sukap.sh", LEAFDIR)<<endl;
        ShellFile<<Form("%s/example/analysis %s %s %d %d %2.2f %s",LEAFDIR,infilename,outfilename,startBin,startBin+neventsPerFile,dr,mPMT?Form("%2.2f",dr_mPMT),"")<<endl;
      }
      ShellFile.close();
    }
  }
}
