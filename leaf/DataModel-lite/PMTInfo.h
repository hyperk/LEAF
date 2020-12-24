#ifndef PMTInfo_H
#define PMTInfo_H

#include <iostream>

#include "Environments.h"
#include "GeoTools.h"

//ROOT Headers
#include "TMath.h"

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
 

class PMTInfo {

	friend class Geometry;

	public:
		PMTInfo(); 
		virtual ~PMTInfo(); 
	 	
	 	int Id_original;
	 	int Id;  
	 	
	 	int Type;

	 	double Position[3]; 
	 	double Orientation[3]; 
	 	
	 	bool mPMT;
	 	
	 	int mPMT_Num;
	 	int mPMT_TubeNum;
	 	int mPMT_Group;
	 	int mPMT_RefTube;
	 	
	 	bool Masked;
	 	
	 	std::vector<double> mPMT_RefX;
	 	std::vector<double> mPMT_RefY;
	 	std::vector<double> mPMT_RefZ;
	 	
	 	int GetTubeNo()		{ return Id;			}
	 	double GetPosition(int i) 	{ return Position[i]; 		} 
	 	double GetOrientation(int i) 	{ return Orientation[i]; 	} 
		
	protected:
		 	
		void Setup_mPMT();
		void Setup_mPMT_Referencial(const PMTInfo lTopPMT);
};

#endif
