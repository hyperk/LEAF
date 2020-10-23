/*********************************************************************************/
/**	HKAstroAnalysis.cc								   **/
/**	Author: Guillaume Pronost (pronost@km.icrr.u-tokyo.ac.jp)		   **/
/**	Date: October 23th 2020						   **/
/**	Desc: Translation of SK Fortran code to C++				   **/
/**	Note: Should not be made open-source without explicit agreement!	   **/
/*********************************************************************************/

#include "HKAstroAnalysis.hh"

HKAstroAnalysis* HKAstroAnalysis::myAna=NULL;


/************************************************************************************************************************/
HKAstroAnalysis::HKAstroAnalysis() {

	myAna = this;
	myManager = HKManager::GetME();
	
	myAna->LoadDirConstant();
	myAna->ResetAnalysis();
}

HKAstroAnalysis::~HKAstroAnalysis() {
	this->ResetAnalysis();
}

HKAstroAnalysis* HKAstroAnalysis::GetME() {

	if ( myAna ) return myAna;

	myAna = new HKAstroAnalysis();
	return myAna;
}
/************************************************************************************************************************/
void HKAstroAnalysis::SetVertex(std::vector<double> lVtx) {

	if ( lVtx.size() == 4 ) {
		fRecoVtx = lVtx;
		
		// Fill variables for additionnal analysis
		fRecoVtx_X	= fRecoVtx[0]; 
		fRecoVtx_Y	= fRecoVtx[1];  
		fRecoVtx_Z	= fRecoVtx[2];  
		fRecoVtx_T	= fRecoVtx[3];  
		fRecoVtx_R2	= fRecoVtx_X * fRecoVtx_X + fRecoVtx_Y * fRecoVtx_Y; 
		fRecoVtx_R	= sqrt(fRecoVtx_R2);
		fRecoVtx_Good	= fRecoVtx[4];  
			
		this->ResetAnalysis();
		
	}
	else {
		std::cout << "HKAstroAnalysis::SetVertex() Invalid vertex set, size= " << lVtx.size() << " should be 4 (X,Y,Z,T)" << std::endl;
	}
}

/************************************************************************************************************************/
void HKAstroAnalysis::SetVertex(float* lVtx) {

	std::vector<double> vVtx;
	vVtx.push_back(lVtx[0]);
	vVtx.push_back(lVtx[1]);
	vVtx.push_back(lVtx[2]);
	vVtx.push_back(lVtx[3]);
	
	this->SetVertex(vVtx);
}

/************************************************************************************************************************/
void HKAstroAnalysis::ResetAnalysis() {

	fInTime20.clear();
	fInTime30.clear();
	fInTime50.clear();
	
	fRecoVtx_Wall 		= 0;
		
	fVtxReco_dir.clear();
	fVtxReco_dirSimple.clear();
	
	fVtxReco_dirKS 	= 0.;
	
	fRecoVtx_E_Simple	= 0.;
	fRecoVtx_E		= 0.;
}
/************************************************************************************************************************/
void HKAstroAnalysis::MakeAnalysis(int iType) {
	this->ResetAnalysis();
	this->ComputeDistanceFromWall();
	
	// Compute In time hit (n50, n30, ...)
	this->ComputeInTimeHit(iType);
			
	if ( fRecoVtx_Wall > 0. && fInTime50.size() > 0 ) {
			
		// Compute simple direction
		// original: 
		//	call lfdir2(bsvertex, bsdir_lfdir2, bsdirks)
		// qmin = .2, qmax = 25., cutoff = .01, twindow = 30.
		fVtxReco_dirSimple = this->GetDirection(0., fInTime30, true, 0.01);
	
		// Compute dirKS from Simple direction
		if ( fVtxReco_dirSimple[3] != -1 ) {
			// Check goodness of direction
			// original:
			//	call lfdirks(bsvertex, bsdir, bsdirks)
			fVtxReco_dirKS = this->GetDirKS(fVtxReco_dirSimple);			
		}				
	}
	else {
		fVtxReco_dirSimple.assign(4,0.);
	}
}

/************************************************************************************************************************/
void HKAstroAnalysis::LoadDirConstant() {

	// From lowe/sklowe/dirlkflf.F
	fEP10[ 0] =  959.1;
	fEP10[ 1] = 1077.29;
	fEP10[ 2] = 1186.07;
	fEP10[ 3] = 1305.84;
	fEP10[ 4] = 1437.70;
	fEP10[ 5] = 1582.88;
	fEP10[ 6] = 1742.72;
	fEP10[ 7] = 2000.  ;
	fEP10[ 8] = 1846.76;
	fEP10[ 9] = 1709.92;
	fEP10[10] = 1584.08;
	fEP10[11] = 1469.25;
	fEP10[12] = 1365.42;
	fEP10[13] = 1272.60;
	fEP10[14] = 1190.78;
	fEP10[15] = 1119.97;
	fEP10[16] = 1060.17;
	fEP10[17] = 1011.37;
	fEP10[18] =  973.57;
	fEP10[19] =  946.78;
	fEP10[20] =  931.  ;
		
	// From lowe/sklowe/lfneweff_sk3_final.F
	fEffHitLambda[ 0] =  7000.;
	fEffHitLambda[ 1] =  8000.;
	fEffHitLambda[ 2] =  9000.;
	fEffHitLambda[ 3] = 10000.;
	fEffHitLambda[ 4] = 11000.;
	fEffHitLambda[ 5] = 12000.;
	fEffHitLambda[ 6] = 13000.;
	fEffHitLambda[ 7] = 14000.;
	fEffHitLambda[ 8] = 15000.;
	fEffHitLambda[ 9] = 16000.;
	fEffHitLambda[10] = 17000.;
	fEffHitLambda[11] = 18000.;
	
	
	// From lowe/sklowe/lfoccor_3.F
	// Occupancy correction table
	// Likely need to be modified to include mPMT later
	int iFullOccupancy  	= 3;
	fNearPMT		= 10;
	for ( int i=0; i < fNearPMT; i++ ) {
		for ( int j=0; j < fNearPMT; j++ ) {
		
			double d1 = (double) i;
			double d2 = (double) j;
			if ( i == j )		
				fOccupancyTable[i][j] = (double) iFullOccupancy;
			else if ( j > i || i == 0 || j == 0 ) 
				fOccupancyTable[i][j] = 0;
			else
				fOccupancyTable[i][j] = - d1 * log(1.0 - d2 / d1) / d2;
				// value: - N/m * log(1 - m/N)
				// 		N: number of total alive PMT
				//		m: number of hit PMT
		}
	}
	
	
}
		
/************************************************************************************************************************/
void HKAstroAnalysis::ComputeInTimeHit(int iType) {
	// Get n hit in a given time window (after tof correction)
	
	//int iNbr = fHitInfo.size();
	int iNbr = myManager->GetNumberOfHit();
	fHitExtInfo.clear();
	
	bool bCond = (iType != AllPMT);
	
	// Compute time of flight and some geometric information
	for ( int iHit=0; iHit < iNbr; iHit++ ) {
		//PMTHit lHit = fHitInfo[iHit];		
		PMTHit lHit = myManager->GetHitInfo(iHit);
		int    iPMT    = lHit.PMT;
		std::vector<double> lPMTInfo = myManager->GetPMTInfo(iPMT);
		
		
		// Skip Hybrid PMT is the option is set
		if ( bCond ) {
			if ( GetPMTType(iPMT) != iType ) continue; 
		}
		PMTHitExt tNewHit;
		
		tNewHit.PMT  = iPMT;
		tNewHit.T    = lHit.T;
		tNewHit.Q    = lHit.Q;
		
		tNewHit.Dist  = GetDistance( lPMTInfo, fRecoVtx );
		tNewHit.ToF   = tNewHit.T - (tNewHit.Dist / CNS2CM);
		
		tNewHit.NormX = (lPMTInfo[0] - fRecoVtx[0])/tNewHit.Dist;
		tNewHit.NormY = (lPMTInfo[1] - fRecoVtx[1])/tNewHit.Dist;
		tNewHit.NormZ = (lPMTInfo[2] - fRecoVtx[2])/tNewHit.Dist;
		
		// Not needed
		tNewHit.theta = 0.;
		tNewHit.phi   = 0.;
		
		// Store hit
		fHitExtInfo.push_back(tNewHit);
	}
	
	iNbr = fHitExtInfo.size();
	
	//std::cout << " InTime computation " << iNbr << " type " << iType << std::endl;
	
	// Clear in time Hit list 
	fInTime20.clear();
	fInTime30.clear();
	fInTime50.clear();

	// Make tmp in time hit list
	std::vector< PMTHitExt > tInTime20;
	std::vector< PMTHitExt > tInTime30;
	std::vector< PMTHitExt > tInTime50;
	
	std::sort(fHitExtInfo.begin(),fHitExtInfo.end(),SortingToF());
	
	for ( int iHit=0; iHit < iNbr; iHit++ ) {
	
		PMTHitExt hHit = fHitExtInfo[iHit];
			
		double dTime = hHit.ToF;
		
		tInTime20 .push_back(hHit);
		tInTime30 .push_back(hHit);
		tInTime50 .push_back(hHit);
		
		while ( (dTime - tInTime20 [0].ToF) > 20. ) {
			tInTime20 .erase(tInTime20 .begin());
		}
		while ( (dTime - tInTime30 [0].ToF) > 30. ) {
			tInTime30 .erase(tInTime30 .begin());
		}
		while ( (dTime - tInTime50 [0].ToF) > 50. ) {
			tInTime50 .erase(tInTime50 .begin());
		}
		
		if ( tInTime20 .size() > fInTime20 .size() )
			fInTime20  = tInTime20;
			
		if ( tInTime30 .size() > fInTime30 .size() )
			fInTime30  = tInTime30;
			
		if ( tInTime50 .size() > fInTime50 .size() )
			fInTime50  = tInTime50;
	}
	
	//std::cout << " InTime computation " << fInTime20.size() << " " << fInTime30.size() << " " << fInTime50.size() << std::endl;
}
/************************************************************************************************************************/
double HKAstroAnalysis::ComputeDistanceFromWall() {
	
	double dR = myManager->GetTankRadius() 	- fRecoVtx_R;
	double dZ = myManager->GetTankHalfHeight() 	- abs(fRecoVtx_Z);
	
	fRecoVtx_Wall = std::min(dR,dZ);
	
	return fRecoVtx_Wall;
}
/************************************************************************************************************************/
double HKAstroAnalysis::GetPMTAngleEfficiency(double dCosPM) {
	// From sklib/coseffsk.F

	if ( dCosPM > 1. ) {
		std::cout << "WARNING: Abnormal dCosPM " << dCosPM << std::endl;
		dCosPM = 1.;	
	}
	if ( dCosPM < 0. ) {
		dCosPM = 0.;	
	}
	
	//double sk1_a[4] = { 0.3542135, 0.5104711, 0.1160429, 0.0118482 }; // SK-I
	double sk2_a[4] = { 0.205349 , 0.523981 , 0.389951 ,-0.131959  }; // SK-II, SK-III, SK-IV, SK-V
	
	double * data = sk2_a;
	
	return (  data[0] 
		+ data[1] * dCosPM
		+ data[2] * dCosPM * dCosPM
		+ data[2] * dCosPM * dCosPM * dCosPM );
}
/************************************************************************************************************************/
std::vector<double> HKAstroAnalysis::TransformInVector(std::vector<double> lDir, double dPhi, double dCosTheta ) {

	std::vector<double> lOut(3,0.);
	
	double dSinTheta = sqrt( 1. - (dCosTheta*dCosTheta) );
	
	double dX = dSinTheta * cos(dPhi);
	double dY = dSinTheta * sin(dPhi);
	double dZ = dCosTheta;
	
	if ( !(abs(lDir[0]) < 1e-2 && abs(lDir[1]) < 1e-2) ) {
	
		double dRDir = sqrt(  	lDir[0] * lDir[0] + 
					lDir[1] * lDir[1] + 
					lDir[2] * lDir[2] );
					
		double dCosA = lDir[2] / dRDir;
		double dSinA = sqrt( 1. - (dCosTheta*dCosTheta) );
		
		double dCosB = lDir[0] / dRDir / dSinA;
		double dSinB = lDir[1] / dRDir / dSinA;
		
		double dXp =  dX * dCosA + dZ * dSinA;
		double dYp =  dY;
		double dZp = -dX * dSinA + dZ * dCosA;
		
		lOut[0] = dXp * dCosB - dYp * dSinB;
		lOut[1] = dXp * dSinB + dYp * dCosB;
		lOut[2] = dZp;
		
		return lOut;
	}
	
	if ( lDir[2] >= 0 ) {
		lOut[0] = dX;
		lOut[1] = dY;
		lOut[2] = dZ;
	}
	else {
		lOut[0] = -dX;
		lOut[1] = -dY;
		lOut[2] = -dZ;
	}
	
	return lOut;
}
/************************************************************************************************************************/
double HKAstroAnalysis::GetDirKS(std::vector<double> lDir) {

	double dCosTh = lDir[2];
	//double dSinTh = sqrt( (1. - dCosTh) * (1. * dCosTh) );
	double dNorm  = sqrt( lDir[0] * lDir[0] + lDir[1] * lDir[1] );
	double dCosPh = lDir[0] / dNorm;
	double dSinPh = lDir[1] / dNorm;
	
	int n50 = fInTime50.size();
	
	std::vector<double> dPhi;
	
	for ( int iHit=0; iHit < n50; iHit++ ) {
	
		double dXX = dSinPh * fInTime50[iHit].NormX - dCosPh * fInTime50[iHit].NormY;
		double dYY = 	  dCosTh * dCosPh * fInTime50[iHit].NormX
				+ dCosTh * dSinPh * fInTime50[iHit].NormY
				-          dSinPh * fInTime50[iHit].NormZ;
				
		dPhi.push_back( atan2( dYY, dXX ) * 180. / 3.1416 + 180. );
	}
	
	std::sort(dPhi.begin(),dPhi.end());
	
	double dNormDeg =  360. / (double) n50;
	double dKSMax   = -360.;
	double dKSMin   =  360.;
	
	for ( int iHit=0; iHit < n50; iHit++ ) {
	
		double dKS = dPhi[iHit] - dNormDeg * (iHit+1);
		
		if ( dKS >= dKSMax ) dKSMax = dKS;
		if ( dKS <= dKSMin ) dKSMin = dKS;
	}
	
	double dDirKS = (dKSMax-dKSMin) / 360.;
	
	return dDirKS;
}
/************************************************************************************************************************/
double HKAstroAnalysis::GetDir_lkflf(double dCos, int iE, bool bSimple) {

	double dOutput = 0.;
	
	if ( bSimple ) {
		// From lowe/sklowe/dirlkflf.F
		
		if ( dCos < -1 || dCos > 1. ) {
			std::cout << "WARNING: dCos out of range in GetDir_lkflf simple " << dCos << std::endl;
			dCos = dCos < 0. ? -1. : 1.; 
		}
		
		if ( dCos <= 0.2 ) {
			dOutput = exp( 5.61 + 1.36 * dCos ) / 2000.;
		}
		else if ( dCos <= 0.6 ) {
			dOutput = exp( 5.42 + 2.41 * dCos ) / 2000.;
		}
		else {
			int    iCos  = (int) ( dCos * 50. ) - 30;
			double dFrac = ( dCos * 50. - 30. ) - (double) iCos;
			
			
			if ( iCos == 20 ) 
				dOutput = fEP10[20];
			else 
				dOutput = fEP10[iCos] + dFrac * ( fEP10[iCos+1] - fEP10[iCos] );
				
			dOutput /= 2000.;
		}
		
		return dOutput;
	}
	
	// From low/sklowe/dirlkflf_ene_sk4.F
	int iCos = (int) ( ( dCos + 1.0 ) * 50 );
	return fHitPat[iE][iCos-1];
}
/************************************************************************************************************************/
std::vector<double> HKAstroAnalysis::GetDirection(double dEnergy, std::vector<PMTHitExt> tInTime, bool bSimple, double dCutOff) {
	// Compute direction and direction prob
	//  Input needed (for dir4): energy, n20, false, 0.2
	//  Input needed (for dir2): energy, n30, true,  0.01
	 
	std::vector<double> tOutput(4,0);
	int nHit = tInTime.size();
	
	int iE = 0;
	
	// Get vector index for GetDir_lkflf
	if ( !bSimple ) {
		if      ( dEnergy <  3.5 ) iE =  0;
		else if ( dEnergy <  4.5 ) iE =  1;
		else if ( dEnergy <  5.5 ) iE =  2;
		else if ( dEnergy <  6.5 ) iE =  3;
		else if ( dEnergy <  7.5 ) iE =  4;
		else if ( dEnergy <  8.5 ) iE =  5;
		else if ( dEnergy <  9.5 ) iE =  6;
		else if ( dEnergy < 10.5 ) iE =  7;
		else if ( dEnergy < 11.5 ) iE =  8;
		else if ( dEnergy < 12.5 ) iE =  9;
		else if ( dEnergy < 13.5 ) iE = 10;
		else if ( dEnergy < 14.5 ) iE = 11;
		else if ( dEnergy < 15.5 ) iE = 12;
		else if ( dEnergy < 16.5 ) iE = 13;
		else if ( dEnergy < 19.0 ) iE = 14;
		else if ( dEnergy < 21.0 ) iE = 15;
		else if ( dEnergy < 23.0 ) iE = 16;
		else if ( dEnergy < 25.0 ) iE = 17;
		else if ( dEnergy < 27.0 ) iE = 18;
		else if ( dEnergy < 40.0 ) iE = 19;
		else if ( dEnergy < 60.0 ) iE = 20;
		else			   iE = 21;
	}
	
	
	if ( nHit > 1 ) {
	
		// Variables:		
		double dDeltaXsum = 0;
		double dDeltaYsum = 0;
		double dDeltaZsum = 0;
	
		for ( int iHit=0; iHit < nHit; iHit++ ) {			
			dDeltaXsum += tInTime[iHit].NormX;
			dDeltaYsum += tInTime[iHit].NormY;
			dDeltaZsum += tInTime[iHit].NormZ;
		}
		
		double dDist = sqrt( 	  dDeltaXsum * dDeltaXsum 
					+ dDeltaYsum * dDeltaYsum 
					+ dDeltaZsum * dDeltaZsum );
					
		tOutput[0] = dDeltaXsum / dDist;
		tOutput[1] = dDeltaYsum / dDist;
		tOutput[2] = dDeltaZsum / dDist;
	
		double dProb = 0;
		
		for ( int iHit=0; iHit < nHit; iHit++ ) {
			int iPMT = tInTime[iHit].PMT;
			std::vector<double> lPMTInfo = myManager->GetPMTInfo(iPMT);
			
			double dCosTh =   tOutput[0] * tInTime[iHit].NormX 
					+ tOutput[1] * tInTime[iHit].NormY 
					+ tOutput[2] * tInTime[iHit].NormZ;
			// 20080718 tentative fix by y.takeuchi  (is this ok??)
			if ( dCosTh > 1. && dCosTh < 1.0001 ) {
				dCosTh = 1.;
			}
			// end
			
			double dTmp = GetDir_lkflf(dCosTh,iE,bSimple);
			if ( dTmp <= dCutOff ) dTmp = dCutOff;
			
			
			double dCosPM = 	- lPMTInfo[3] * tInTime[iHit].NormX
						- lPMTInfo[4] * tInTime[iHit].NormY
						- lPMTInfo[5] * tInTime[iHit].NormZ;
					
			//double dRR    = GetDistance(lPMTInfo,fRecoVtx);
			
			// 20080718 tentative fix by y.takeuchi  (is this ok??)
			if ( dCosPM > 1. && dCosPM < 1.0001 ) {
				dCosPM = 1.;
			}
			// end
			
			double dCorEff = dCosPM / GetPMTAngleEfficiency(dCosPM) * tInTime[iHit].Q;
			
			dProb += (log10(dTmp) * dCorEff);
		}
		// Save prob
		tOutput[3] = dProb;
		
		double dTheta[4] = 	{
						20. *3.1416/180.,
						 9. *3.1416/180.,
						 4. *3.1416/180.,
						 1.6*3.1416/180.,
					};
					
		for ( int iStep = 0; iStep < 4; iStep++ ) {
			
			double dCosTheta = cos(dTheta[iStep]);
			
			int iTry  = 0;
			int iMove = 1;
			
			while ( iTry < 9 && iMove == 1 ) {
				iMove  = 0;
				iTry  += 1;
				
				for ( int iPhi=0; iPhi < 6; iPhi++ ) {
					
					double dPhi = float(iPhi) / 6. * 2. * 3.1416;
					
					std::vector<double> lWDir = TransformInVector(tOutput,dPhi,dCosTheta);
					
					double dProb2 = 0.;
					
					for ( int iHit=0; iHit < nHit; iHit++ ) {
						int iPMT = tInTime[iHit].PMT;
						std::vector<double> lPMTInfo = myManager->GetPMTInfo(iPMT);
						
						double dCosTh =   lWDir[0] * tInTime[iHit].NormX
								+ lWDir[1] * tInTime[iHit].NormY 
								+ lWDir[2] * tInTime[iHit].NormZ;
								
						// 20080718 tentative fix by y.takeuchi  (is this ok??)
						if ( dCosTh > 1. && dCosTh < 1.0001 ) {
							dCosTh = 1.;
						}
						// end
						
						double dTmp = GetDir_lkflf(dCosTh,iE,bSimple);
						if ( dTmp <= dCutOff ) dTmp = dCutOff;
												
						double dCosPM = 	- lPMTInfo[3] * tInTime[iHit].NormX
									- lPMTInfo[4] * tInTime[iHit].NormY
									- lPMTInfo[5] * tInTime[iHit].NormZ;
								
						//double dRR    = GetDistance(lPMTInfo,fRecoVtx);
						
						// 20080718 tentative fix by y.takeuchi  (is this ok??)
						if ( dCosPM > 1. && dCosPM < 1.0001 ) {
							dCosPM = 1.;
						}
						// end			
						
						double dCorEff = dCosPM / GetPMTAngleEfficiency(dCosPM) * tInTime[iHit].Q;
						
						dProb2 += (log10(dTmp) * dCorEff);
					}
					
					if ( dProb2 > dProb ) {
						dProb   = dProb2;
						iMove   = 1;
						tOutput = lWDir;
						tOutput.push_back( dProb );
					}	
				} // for iPhi
			} // while
		} // for iStep
	}
	
	return tOutput;
}
/************************************************************************************************************************/
double HKAstroAnalysis::GoodnessBonsai() {
	// From timefit.cc (bonsai code)
	// fRecoVtx
	
	double g, count;
	
	g     = 0.;
	count = 0.;
	
	for ( unsigned int i=0; i < fHitExtInfo.size(); i++ ) {
		double dDeltaT = fHitExtInfo[i].ToF - fRecoVtx[3];
		
		double dt      = dDeltaT * dDeltaT * 0.04;
		double w       = dt * 0.0035;
		
		//std::cout << " Hit " << i << " ToF - VtxTime " << dDeltaT << " " << w << " " << dt << std::endl;
		if ( w > 18. ) continue;
		
		w = exp(-w);
		count += w;
		dt *= 0.5;
		
		if ( dt <= 50 ) 
			g += w * exp(-dt);
	}
	
	//std::cout << " g/count " << g << " / " << count << " = " << g/count << std::endl;
	if ( count == 0 ) {
		return 0;
	}
	
	return g/count;
}
