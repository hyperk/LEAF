#ifndef RECOINFO_H
#define RECOINFO_H

#include <iostream>
#include <vector>
#include <cmath>

#include "Environments.h"
#include "TimeDelta.h"

#include "RootHitInfo.h"


class RootRecoInfo {

	public:
	
		RootRecoInfo();
		virtual ~RootRecoInfo();
		
		// Vertex fitter algorithm 
		// Base of the RecoInfo
		HKAA::RecoVtxInfo AlgVtx; 
		
		std::vector<double> 			Vtx;
		Pos3D					Pos;
		TimeDelta				Time;
			
		double					Goodness;	
		double					Goodness_Time;
		
		double					Wall; // Distance from wall
	
		std::vector<double> 			E;
		std::vector< std::vector<double> > 	Dir;
		std::vector< std::vector<double> > 	Cone;
		std::vector<double>			Dir_Goodness;
		
		// For AstroAnalysis:
		RootHitCollection			ExtendedHitCollection;
		RootHitCollection			InTime20;
		RootHitCollection			InTime30;
		RootHitCollection			InTime50;
		
		//
			  
		void Reset();

	private:

		ClassDef(RootRecoInfo,1) //EventRootInfo structure

};

#if !defined(__CLING__)
ClassImp(RootRecoInfo)
#endif

#endif //RECOINFO_H
