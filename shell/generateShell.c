{
  char * ShellName = new char[256];
  ofstream ShellFile(ShellName);
  char * infilename = new char[512];
  char * outfilename = new char[512];
  char * LEAFDIR = getenv("LEAFDIR");
  
  int nevents=10000;
  int neventsPerFile=20;
  for(int nrj=10;nrj<11;nrj++){
  //for(int nrj=3;nrj<16;nrj++){
    sprintf(infilename,"%s/inputs/wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_fulltank_0hitstrigger_10000.root",LEAFDIR,nrj);
    
    for(int startBin=0;startBin<nevents;startBin+=neventsPerFile){
      sprintf(ShellName,"%s/jobs/fitLE_e%d_start%d.sh",LEAFDIR,nrj,startBin);
      sprintf(outfilename,"%s/outputs/LEAF_test_%dMeV_tmp%d.root",LEAFDIR,nrj,startBin);

      ofstream ShellFile(ShellName);
      if(ShellFile){
	ShellFile<<"#!/bin/bash"<<endl;
	ShellFile<<Form("%s/app/FitVertexLE -f %s -o %s -s %d -e %d",LEAFDIR,infilename,outfilename,startBin,startBin+neventsPerFile)<<endl;
      }
    }
  }
}
