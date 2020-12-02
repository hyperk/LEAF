#ifndef PMTInfo_H
#define PMTInfo_H

#include <iostream>

#include "Environments.h"
#include "GeoTools.h"

//ROOT Headers
#include "TMath.h"
#include "TObject.h"

/**
* \class PMTInfo
 *
 * This class holds PMT geometry information (PMT type, position, direction, etc.)
 *
 *
 * $Author: G.Pronost $
 * $Date: 2020/11/13 16:15:00 $
 * Contact: pronost@km.icrr.u-tokyo.ac.jp
 *          
 */
 

class RootPMTInfo : public TObject {

	friend class Geometry;

	public:
		RootPMTInfo(); 
	 	
	 	int Id_original;
	 	int Id;  
	 	
	 	int Type;

	 	double Position[3]; 
	 	double Orientation[3]; 
	 	
	 	bool mPMT;
	 	
	 	int mPMT_TubeNum;
	 	int mPMT_Group;
	 	int mPMT_RefTube;
	 	
	 	std::vector<double> mPMT_RefX;
	 	std::vector<double> mPMT_RefY;
	 	std::vector<double> mPMT_RefZ;
		
	private:
		 	
		void Setup_mPMT();
		void Setup_mPMT_Referencial(const RootPMTInfo lTopPMT);
		
		ClassDef(RootPMTInfo,1) //EventRootInfo structure
};

#if !defined(__CLING__)
ClassImp(RootPMTInfo)
#endif


#endif
