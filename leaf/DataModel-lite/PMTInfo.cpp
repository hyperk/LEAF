#include "PMTInfo.h"

PMTInfo::PMTInfo() {
	 
	Id_original	= 0;
	Id		= 0;
	Type		= HKAA::kIDPMT_BnL;
	
	// PMT Position initialization
	Position[0]	= 0.;
	Position[1]	= 0.;
	Position[2]	= 0.;
	
	// PMT Direction initialization
	Orientation[0]	= 0.;
	Orientation[1]	= 0.;
	Orientation[2]	= 0.;
	
	// mPMT related info initialization
	mPMT 		= false;
	mPMT_TubeNum 	= 0;
	mPMT_Group 	= 0;
	mPMT_RefTube 	= 0;
	
	// Mask information
	Masked 	= true;
}

PMTInfo::~PMTInfo() {
	mPMT_RefX.clear();
	mPMT_RefY.clear();
	mPMT_RefZ.clear();
}

void PMTInfo::Setup_mPMT() {

	mPMT_RefTube = Id + ( HKAA::kmPMT_TopID - mPMT_TubeNum );
	
	if ( mPMT_TubeNum <= 12 )	mPMT_Group = 0;
	else if ( mPMT_TubeNum <= 18)	mPMT_Group = 1;
	else				mPMT_Group = HKAA::kmPMT_Groups - 1;
}

void PMTInfo::Setup_mPMT_Referencial(const PMTInfo lTopPMT) {
	// Function from LEAF
	
	if ( lTopPMT.Id != mPMT_RefTube ) {
		std::cout << " Incorrect reference tube send to setup mPMT Referencial " << std::endl; 
		return;
	}
	
	mPMT_RefX.resize(3.,0);
	mPMT_RefY.resize(3.,0);
	mPMT_RefZ.resize(3.,0);
	
	double dEX[3];
	double dEY[3];
	double dEZ[3];
	
	//1. Get the PMT Orientation
	dEZ[0] = Orientation[0];
	dEZ[1] = Orientation[1];
	dEZ[2] = Orientation[2];
	
	// Save in global
  	mPMT_RefZ[0] = dEZ[0];
  	mPMT_RefZ[1] = dEZ[1];
  	mPMT_RefZ[2] = dEZ[2];
  	
	/*
	if ( VERBOSE >= 3 ) {
		std::cout << " PMT number = " << Id << ", Tube number in mPMT = " << mPMT_TubeNum << std::endl; 
		std::cout << " Position of the PMT = " << Position[0] << ", " << Position[1] << ", " << Position[1] << std::endl;
		std::cout << " Orientation of the PMT = " << mPMT_RefZ[0] << ", " << mPMT_RefZ[1] << ", " << mPMT_RefZ[2] << std::endl;
	}
	*/
	
	//2. Get the mPMT position and its direction -> Provides orgin and first axis of the referential (ez)
	//We do not have the mPMT position, so we will take the position of the PMT on top of the mPMT: number:
		
	/*
	if ( VERBOSE >= 3 ) {
		std::cout << " Top PMT number = " << lTopPMT.Id << ", Tube number in mPMT = " << lTopPMT.mPMT_TubeNum << std::endl; 
		std::cout << " Position of the PMT = " << lTopPMT.Position[0] << ", " << lTopPMT.Position[1] << ", " << lTopPMT.Position[2] << std::endl;
		std::cout << " Orientation of the PMT = " << lTopPMT.Orientation[0] << ", " << lTopPMT.Orientation[1] << ", " << lTopPMT.Orientation[2] << std::endl;
	}
	*/	

	//3. Define the second vector of the referential in the orthogonal plane to ez and colinear to the center PMT -> Active PMT plane
	//a. PMT position in the mPMT referencec plane
	
	double dPosition_Ref[3];
	
	dPosition_Ref[0] = Position[0] - lTopPMT.Position[0];
	dPosition_Ref[1] = Position[1] - lTopPMT.Position[1];
	dPosition_Ref[2] = Position[2] - lTopPMT.Position[2];
	
	GeoTools::Normalize(dPosition_Ref);
	
	/*
	if ( VERBOSE >= 3 ) {
		std::cout << " PMT position in the mPMT referential = " << dPosition_Ref[0] << ", " << dPosition_Ref[1] << ", " << dPosition_Ref[2] << std::endl;
	}
	*/
	
	//b. Define ey, perpendicular to ez and the center PMT -> Active PMT vector  
	
	GeoTools::CrossProduct(dEZ,dPosition_Ref,dEY);
	GeoTools::Normalize(dEY);
	
	// Save in global
  	mPMT_RefY[0] = dEY[0];
  	mPMT_RefY[1] = dEY[1];
  	mPMT_RefY[2] = dEY[2];
  	/*
	if ( VERBOSE >= 3 ) {
		std::cout << " Scal product test = " 	<< TMath::ACos(Astro_GetScalarProd(mPMT_RefY,mPMT_RefZ))*180/TMath::Pi() << ", " 
							<< Astro_GetLength(mPMT_RefY) << ", " << Astro_GetLength(mPMT_RefZ) << std::endl;
	}
	*/
  	
	//4. And then create ex which should be orthogonal to the 2 remaining vectors
	GeoTools::CrossProduct(dEY,dEZ,dEX);
	
	// Save in global
  	mPMT_RefX[0] = dEX[0];
  	mPMT_RefX[1] = dEX[1];
  	mPMT_RefX[2] = dEX[2];
  	
	//As ez and ey should be unitary, so does ex. So, let's check as a debug
  	
  	double dLengthX = Astro_GetLength(mPMT_RefX);
  	// dLengthX is not exactly 1 let's allow +/- 1e-6
  	
  	if ( dLengthX > (1 + 1e-6) || dLengthX < (1 - 1e-6) ) {
  		std::cout << " There is a referential vector unitarity issue for PMT " << Id << ", length is = "<< dLengthX << std::endl;
  		
		std::cout << " \t Top PMT number = " << lTopPMT.Id << ", Tube number in mPMT = " << lTopPMT.mPMT_TubeNum << std::endl; 
		std::cout << " \t PMT position in the mPMT referential = " << dPosition_Ref[0] << ", " << dPosition_Ref[1] << ", " << dPosition_Ref[2] << std::endl;
		std::cout << " \t Scalar product test = " 	<< TMath::ACos(Astro_GetScalarProd(mPMT_RefY,mPMT_RefZ))*180/TMath::Pi() << ", Y: " 
								<< Astro_GetLength(mPMT_RefY) << ", Z: " << Astro_GetLength(mPMT_RefZ) << std::endl;
  	}
  	
	//Conclusion: we now have our referential.
}
