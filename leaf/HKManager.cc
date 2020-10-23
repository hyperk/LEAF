/*********************************************************************************/
/**	HKManager.cc							   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: June 17th 2020							   **/
/**	Desc: Manage fitter and PMT hits for HK analysis			   **/
/*********************************************************************************/

#include "HKManager.hh"

HKManager* HKManager::myManager=NULL;

/************************************************************************************************************************/

HKManager::HKManager() {
	myManager = this;
}

HKManager::~HKManager() {

	fPMT_Info	.clear();
	fPMT_Group	.clear();
	fPMT_TubeInMPMT	.clear();
	fPMT_RefX	.clear();
	fPMT_RefY	.clear();
	fPMT_RefZ	.clear();
	fPMT_RefInMPMT	.clear();
	fHitInfo	.clear();
	
	myManager = NULL;
	
	
	fPMTDiagDistance  = 0;
	fPMTBarrelExtremZ_Lo = 0;
	fPMTBarrelExtremZ_Hi = 0;
	
	fNumberOfPMT = 0;
	fNumberOfmPMT = 0;
}

HKManager* HKManager::GetME() {

	if ( myManager ) return myManager;

	myManager = new HKManager();
	return myManager;
}

/************************************************************************************************************************/

void HKManager::SetGeometry( WCSimRootGeom * wGeo, double dDarkRate_Normal, double dDarkRate_mPMT ) {

	// Keep RootGeometry event:
	fWCGeo = wGeo;
	
	// Read important variables:
	fTankRadius = fWCGeo->GetWCCylRadius(); //Conversion from cm to m removed
	fTankHeight = fWCGeo->GetWCCylLength(); //Conversion from cm to m removed
	fTankHalfHeight = fTankHeight / 2.;
	
	fDarkRate_Normal = dDarkRate_Normal;
	fDarkRate_mPMT   = dDarkRate_mPMT;

	if ( dDarkRate_Normal == 0 ) {
#ifdef WCSIM_single_PMT_type	
		fDarkRate_ns[0] = 4200.  * 2e4 * /*fWCGeo->GetWCNumPMT() */ 1e-9;
		fDarkRate_ns[1] = 0;
#else
		fDarkRate_ns[0] = 8400.  * 2e4 * /*fWCGeo->GetWCNumPMT(false) */ 1e-9;
		fDarkRate_ns[1] = 100.   * 5e3 * /*fWCGeo->GetWCNumPMT(true ) */ 1e-9 * mPMT_TOP_ID;
#endif
	}
	else {	
#ifdef WCSIM_single_PMT_type	
		fDarkRate_ns[0] = fDarkRate_Normal * fWCGeo->GetWCNumPMT() * 1e-9;
		fDarkRate_ns[1] = 0;
#else

		int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
		int iNbr_mPMT = fWCGeo->GetWCNumPMT(true );
#ifdef BUG_WCGEO
		if ( iNbr_Norm > 40e3 ) {
			iNbr_mPMT = iNbr_Norm;
		}
		else {
			// Likely no mPMT (or maybe correct output)
			iNbr_mPMT = fWCGeo->GetWCNumPMT(true);
		}
#endif	

		fDarkRate_ns[0] = fDarkRate_Normal * iNbr_Norm * 1e-9;
		fDarkRate_ns[1] = fDarkRate_mPMT   * iNbr_mPMT * 1e-9 * mPMT_TOP_ID;
#endif
	}
	
	this->LoadPMTInfo();
}

/************************************************************************************************************************/
// Geometry functions:

void HKGeoTools::CrossProduct( double* a, double* b, double* c ) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
} 

void HKGeoTools::Normalize( double* a ) {

	double dL = GetLength(a);
	
	a[0] /= dL;
	a[1] /= dL;
	a[2] /= dL;
} 

double HKManager::GroupPMTs(int pmtType, int pmt_number) {
	if(pmtType==0) 			return 0;
	else{
		if(pmt_number<=12) 	return 0;
		else if(pmt_number<=18) return 1;
		else 			return NGROUP_PMT - 1;
	}
}

void HKManager::LoadPMTInfo() {

	// clear
	fPMT_TubeInMPMT.clear();
	fPMT_Group.clear();
	fPMT_Info.clear();
	
#ifdef WCSIM_single_PMT_type	
	int iNbr_Norm = fWCGeo->GetWCNumPMT();
	int iNbr_mPMT = 0;
#else
	int iNbr_Norm = fWCGeo->GetWCNumPMT(false);
	int iNbr_mPMT = fWCGeo->GetWCNumPMT(true);
	
	
#ifdef BUG_WCGEO
	if ( iNbr_mPMT == iNbr_Norm ) {
		iNbr_mPMT = 0;
	}
	
	if ( iNbr_Norm > 40e3 ) {
		iNbr_mPMT = iNbr_Norm;
	}
	else {
		// Likely no mPMT (or maybe correct output)
	}
	
#endif	
#endif
	std::cout << " LEAF setting # PMT = " << iNbr_Norm << std::endl;
	std::cout << " LEAF setting # mPMT = " << iNbr_mPMT << std::endl;
	
	fPMT_Info	.resize(MAX_PMT);
	fPMT_Group	.assign(MAX_PMT,0);
	fPMT_TubeInMPMT	.assign(MAX_PMT,0);
	fPMT_RefX	.resize(MAX_PMT);
	fPMT_RefY	.resize(MAX_PMT);
	fPMT_RefZ	.resize(MAX_PMT);
	fPMT_RefInMPMT	.assign(MAX_PMT,0);
	
	
	WCSimRootPMT wPMT;
	
	double dMin1 = 0;
	double dMin2 = 0;
	
	fNumberOfPMT = 0;
	// Normal PMTs
	for ( int iPMT=0; iPMT < iNbr_Norm; iPMT++ ) {
	
#ifdef WCSIM_single_PMT_type	
		wPMT = fWCGeo->GetPMT(iPMT);
#else
		wPMT = fWCGeo->GetPMT(iPMT,false);
#endif
		std::vector<double> vPMT(6,0);
		
		
		for ( int j=0; j < 3; j++ ) {
			vPMT[j  ] = wPMT.GetPosition(j);
			vPMT[j+3] = wPMT.GetOrientation(j);
		}
		
		if ( vPMT[2] < dMin1 ) {
			dMin2 = dMin1;
			dMin1 = vPMT[2];
		}
		
		if ( vPMT[2] > dMin1 && vPMT[2] < dMin2 ) {
			dMin2 = vPMT[2];
		}
		
		fPMT_TubeInMPMT[iPMT] = 0;
		fPMT_Group     [iPMT] = 0;
		fPMT_Info      [iPMT] = vPMT;
		
		fPMT_RefInMPMT [iPMT] = 0;
		
		if ( !( vPMT[0] == 0 && vPMT[1] == 0 && vPMT[2] == 0 ) ) {	
			fNumberOfPMT += 1;
		}
	}
	
	
	// Compute minimal diagonal distance:
	// PMT 4 is the closest diag PMT from PMT 0
	
	fPMTDiagDistance = std::ceil( sqrt( pow(fPMT_Info[0][0]-fPMT_Info[4][0],2.) + pow(fPMT_Info[0][1]-fPMT_Info[4][1],2.) + pow(fPMT_Info[0][2]-fPMT_Info[4][2],2.) ) );
	fPMTDiagDistance += 2.; // Add 2cm for safety
	
	fPMTBarrelExtremZ_Lo = std::abs( std::ceil(dMin2) ) + 1;
	fPMTBarrelExtremZ_Hi = std::abs( std::ceil(dMin1) ) + 1;
	
	
	
	
#ifndef WCSIM_single_PMT_type	
	// mPMTs
	fNumberOfmPMT = 0;
	for ( int iPMT=0; iPMT < iNbr_mPMT; iPMT++ ) {
	
		wPMT = fWCGeo->GetPMT(iPMT,true);
		
		std::vector<double> vPMT(6,0);
		
		for ( int j=0; j < 3; j++ ) {
			vPMT[j  ] = wPMT.GetPosition(j);
			vPMT[j+3] = wPMT.GetOrientation(j);
		}
		
		int iNewID = iPMT+mPMT_ID_SHIFT;
		
		int iPMTNo = wPMT.GetmPMT_PMTNo();
		
		if ( iPMTNo < 0 || iPMTNo > 30 ) {
			if ( iPMT > 0 ) {
				std::cout << " ERROR mPMT " << iPMT << "("<< iNewID << ") is not defined (" << vPMT[0] << ", " << vPMT[1] << ", " << vPMT[2] 
					<< ") (" << vPMT[3] << ", " << vPMT[4] << ", " << vPMT[5] << ") " << iPMTNo << std::endl;
			}
			iPMTNo = 0;
			
		}
		
		fPMT_TubeInMPMT[iNewID] = iPMTNo;
		fPMT_Group     [iNewID] = this->GroupPMTs(1,iPMTNo);
		fPMT_Info      [iNewID] = vPMT;
		
		fPMT_RefInMPMT [iNewID] = iNewID + ( mPMT_TOP_ID - iPMTNo );
		
		if ( fPMT_RefInMPMT[iNewID] >= (iNbr_mPMT+mPMT_ID_SHIFT) ) {
			std::cout << " ERROR: Reference PMT for mPMT is outside mPMT range " << std::endl;
		}
		
		if ( !( vPMT[0] == 0 && vPMT[1] == 0 && vPMT[2] == 0 ) ) {	
			fNumberOfmPMT += 1;
		}
	}
	
	// Re-loop to make PMT Referencial
	for ( int iPMT=0; iPMT < fNumberOfmPMT; iPMT++ ) {
		int iNewID = iPMT+mPMT_ID_SHIFT;
		this->MakeMPMTReferencial(iNewID);
	}
	
#endif
}

void HKManager::MakeMPMTReferencial(int iPMT) {

	if ( !GetPMTType(iPMT) ) return;

	// fPMT_RefX[iPMT] 
	// fPMT_RefY[iPMT]
	// fPMT_RefZ[iPMT]
	
	fPMT_RefX[iPMT].resize(3.,0);
	fPMT_RefY[iPMT].resize(3.,0);
	fPMT_RefZ[iPMT].resize(3.,0);
	
	double dEX[3];
	double dEY[3];
	double dEZ[3];
	
	int iTubeInMPMT = fPMT_TubeInMPMT[iPMT];
  
	//1. Get the PMT position
	// -> fPMT_Info[iPMT][0]
	// -> fPMT_Info[iPMT][1]
	// -> fPMT_Info[iPMT][2]
	// -> fPMT_Info[iPMT][3] = ez[0]
	// -> fPMT_Info[iPMT][4] = ez[1]
	// -> fPMT_Info[iPMT][5] = ez[2]
	
	dEZ[0] = fPMT_Info[iPMT][3];
	dEZ[1] = fPMT_Info[iPMT][4];
	dEZ[2] = fPMT_Info[iPMT][5];
	
	// Save in global
  	fPMT_RefZ[iPMT][0] = dEZ[0];
  	fPMT_RefZ[iPMT][1] = dEZ[1];
  	fPMT_RefZ[iPMT][2] = dEZ[2];
	
	if ( VERBOSE >= 3 ) {
		std::cout << " PMT number = " << iPMT << ", Tube number in mPMT = " << iTubeInMPMT << std::endl; 
		std::cout << " Position of the PMT = " << fPMT_Info[iPMT][0] << ", " << fPMT_Info[iPMT][1] << ", " << fPMT_Info[iPMT][1] << std::endl;
		std::cout << " Orientation of the PMT = " << fPMT_RefZ[iPMT][0] << ", " << fPMT_RefZ[iPMT][1] << ", " << fPMT_RefZ[iPMT][2] << std::endl;
	}
	
	//2. Get the mPMT position and its direction -> Provides orgin and first axis of the referential (ez)
	//We do not have the mPMT position, so we will take the position of the PMT on top of the mPMT: number:
	
	int iPMTTop = fPMT_RefInMPMT[iPMT];
	
	if ( VERBOSE >= 3 ) {
		std::cout << " Top PMT number = " << iPMTTop << ", Tube number in mPMT = " << fPMT_TubeInMPMT[iPMTTop] << std::endl; 
		std::cout << " Position of the PMT = " << fPMT_Info[iPMTTop][0] << ", " << fPMT_Info[iPMTTop][1] << ", " << fPMT_Info[iPMTTop][2] << std::endl;
		std::cout << " Orientation of the PMT = " << fPMT_Info[iPMTTop][3] << ", " << fPMT_Info[iPMTTop][4] << ", " << fPMT_Info[iPMTTop][5] << std::endl;
	}
	

	//3. Define the second vector of the referential in the orthogonal plane to ez and colinear to the center PMT -> Active PMT plane
	//a. PMT position in the mPMT referencec plane
	
	double dPosition_Ref[3];
	
	dPosition_Ref[0] = fPMT_Info[iPMT][0] - fPMT_Info[iPMTTop][0];
	dPosition_Ref[1] = fPMT_Info[iPMT][1] - fPMT_Info[iPMTTop][1];
	dPosition_Ref[2] = fPMT_Info[iPMT][1] - fPMT_Info[iPMTTop][2];
	
	HKGeoTools::Normalize(dPosition_Ref);
	
	if ( VERBOSE >= 3 ) {
		std::cout << " PMT position in the mPMT referential = " << dPosition_Ref[0] << ", " << dPosition_Ref[1] << ", " << dPosition_Ref[2] << std::endl;
	}
	
	//b. Define ey, perpendicular to ez and the center PMT -> Active PMT vector  
	
	HKGeoTools::CrossProduct(dEZ,dPosition_Ref,dEY);
	HKGeoTools::Normalize(dEY);
	
	// Save in global
  	fPMT_RefY[iPMT][0] = dEY[0];
  	fPMT_RefY[iPMT][1] = dEY[1];
  	fPMT_RefY[iPMT][2] = dEY[2];
  	
	if ( VERBOSE >= 3 ) {
		std::cout << " Scal product test = " 	<< TMath::ACos(GetScalarProd(fPMT_RefY[iPMT],fPMT_RefZ[iPMT]))*180/TMath::Pi() << ", " 
							<< GetLength(fPMT_RefY[iPMT]) << ", " << GetLength(fPMT_RefZ[iPMT]) << std::endl;
	}
  	
	//4. And then create ex which should be orthogonal to the 2 remaining vectors
	HKGeoTools::CrossProduct(dEY,dEZ,dEX);
	
	// Save in global
  	fPMT_RefX[iPMT][0] = dEX[0];
  	fPMT_RefX[iPMT][1] = dEX[1];
  	fPMT_RefX[iPMT][2] = dEX[2];
  	
	//As ez and ey should be unitary, so does ex. So, let's check as a debug
  	
  	double dLengthX = GetLength(fPMT_RefX[iPMT]);
  	// dLengthX is not exactly 1 let's allow +/- 1e-6
  	if ( dLengthX > (1 + 1e-6) || dLengthX < (1 - 1e-6) ) {
  		std::cout << " There is a referential vector unitarity issue for PMT " << iPMT << ", length is = "<< dLengthX << std::endl;
  		
		std::cout << " \t Top PMT number = " << iPMTTop << ", Tube number in mPMT = " << fPMT_TubeInMPMT[iPMTTop] << std::endl; 
		std::cout << " \t PMT position in the mPMT referential = " << dPosition_Ref[0] << ", " << dPosition_Ref[1] << ", " << dPosition_Ref[2] << std::endl;
		std::cout << " \t Scalar product test = " 	<< TMath::ACos(GetScalarProd(fPMT_RefY[iPMT],fPMT_RefZ[iPMT]))*180/TMath::Pi() << ", Y: " 
								<< GetLength(fPMT_RefY[iPMT]) << ", Z: " << GetLength(fPMT_RefZ[iPMT]) << std::endl;
  	}
  	
  	
	//Conclusion: we now have our referential.

}
