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
  //for(int nrj=10;nrj<11;nrj++){
  for(int nrj=3;nrj<16;nrj++){
    //sprintf(infilename,"/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_closewall_0hitstrigger_100000.root",nrj);
    sprintf(infilename,"/disk01/usr5/bquilain/wcsim_hkhybridmpmt10pc14374100Hz_4.2kHzBL_e%d_center_nominal_fulltank_0hitstrigger_100000.root",nrj);

    for(int startBin=0;startBin<nevents;startBin+=neventsPerFile){
      sprintf(ShellName,"%s/jobs/fitLE_e%d_start%d.sh",LEAFDIR,nrj,startBin);
      sprintf(outfilename,"/disk01/usr5/bquilain/LEAFv2_hkhybridmpmt10pc14374100Hz_4.2kHzbl_%dMeV_tmp%d.root",nrj,startBin);
    
      ofstream ShellFile(ShellName);
      if(ShellFile){
	ShellFile<<"#!/bin/bash"<<endl;
	ShellFile<<Form("rm -f %s",outfilename)<<endl;
	ShellFile<<Form("source /disk01/usr5/bquilain/LEAF/LEAF_v1.5/LEAF/RunAtStart_sukap.sh")<<endl;
	ShellFile<<Form("%s/app/analysis %s %s %d %d %2.2f %s",LEAFDIR,infilename,outfilename,startBin,startBin+neventsPerFile,dr,mPMT?Form("%2.2f",dr_mPMT),"")<<endl;
	//ShellFile<<Form("/home/bquilain/nuPRISM/BenjaminAnalysis/FitVertexLE_v2 -f %s -o %s -s %d -e %d",infilename,outfilename,startBin,startBin+neventsPerFile)<<endl;
	//ShellFile<<Form("%s/app/FitVertexLE -f %s -o %s -s %d -e %d",LEAFDIR,infilename2,outfilename2,startBin,startBin+neventsPerFile)<<endl;
      }
    }
  }
}
