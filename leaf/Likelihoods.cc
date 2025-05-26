#include "Likelihoods.hh"
#include "LeafInputs.hh"

//Function to optimize, for MIGRAD
void MinuitLikelihood(int & /*nDim*/, double * /*gout*/, double &NLL, double par[], int /*flg*/)
{
	std::vector<double> vertexPosition(4, 0.); // In centimeters
	for (int i = 0; i < 4; i++) vertexPosition[i] = par[i];
	int nhits = par[5];
	double lowerLimit = par[6];
	double upperLimit = par[7];
	double directionality = par[9];
	// std::vector<double> vertexDirection(3,0.);

	NLL = Likelihoods::Vertex_Time_NLL(fHitCollection, vertexPosition, nhits, lowerLimit, upperLimit, true, false, directionality);
}

//Function to optimize, for MIGRAD
void MinuitDirNLL(int& /*nDim*/, double* /*gout*/, double& DNLL, double par[], int /*flg*/) 
{
	// Extract theta and phi from parameters
	double theta = par[0];
	double phi = par[1];
	int nhits = static_cast<int>(par[2]);

	// Ensure the vertexPosition vector is properly initialized
	std::vector<double> vertexPosition(4, 0.0);
	for (int i = 0; i < 4; i++) vertexPosition[i] = par[i + 3];

	// Ensure theta and phi are within valid ranges
	if (theta < 0 || theta > TMath::Pi()) 
	{
		DNLL = 1e10; // Assign a large NLL to invalid directions
		return;
	}
	if (phi < -TMath::Pi() || phi > TMath::Pi()) 
	{
		DNLL = 1e10;
		return;
	}

	// Calculate NLL using theta and phi
	DNLL = Likelihoods::Dir_NLL(fHitCollection, vertexPosition, theta, phi, nhits);
}

void MinuitJointNLL(int& nDim, double * gout, double & NLL, double par[], int flg)
{
	std::vector<double> vertexPosition(4, 0.);
	for (int i = 0; i < 4; i++) vertexPosition[i] = par[i];
	double theta = par[4];
	double phi = par[5];
	int nhits = par[7];
	double lowerLimit = par[8];
	double upperLimit = par[9];
	// double directionality = par[9];

	// Ensure theta and phi are within valid ranges
	if (theta < 0 || theta > TMath::Pi()) 
	{
		NLL = 1e10; // Assign a large NLL to invalid directions
		return;
	}
	if (phi < -TMath::Pi() || phi > TMath::Pi()) 
	{
		NLL = 1e10;
		return;
	}

	double VtxNLL = Likelihoods::Vertex_Time_NLL(fHitCollection, vertexPosition, nhits, lowerLimit, upperLimit, true, false, 0);
	double DirNLL = Likelihoods::Dir_NLL(fHitCollection, vertexPosition, theta, phi, nhits);

	NLL = VtxNLL * DirNLL;
}

double Likelihoods::Vertex_Time_NLL(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality)
{
	double NLL = 0;

	for (int ihit = 0; ihit < nhits; ihit++)
	{
		Hit lHit = lHitCol->At(ihit);

		int iPMT = lHit.PMT;
		// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];

		int pmtType = Astro_GetPMTType(iPMT);
		// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

		// double tof = distance / fLightSpeed;
		// double residual = hitTime - tof - vertexPosition[3];
		double residual = ComputeResidualTime(vertexPosition, vertexPosition[3], lHitCol->At(ihit));

		residual = ComputeResidualTime(vertexPosition, vertexPosition[3], lHit);

		double proba = 0;

		bool condition = residual > lowerLimit && residual < upperLimit;

		if (condition)
		{
			proba = fSplineTimePDFQueue[pmtType]->Eval(residual);

#ifdef VERBOSE_NLL
			if (VERBOSE >= 3)
			{
				std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << std::endl;
				std::cout << "Residual=" << residual << ", proba=" << proba << ", pmt type=" << pmtType << std::endl;
			}
#endif

			if (scaleDR)
			{
				std::cout << "scaleDR with Lower Limit: " << lowerLimit << " Upper Limit: " << upperLimit << std::endl;
				EventInfo fEventInfo = MakeEventInfo(lowerLimit, upperLimit, pmtType);

				// First, substract the DR from the PDF to keep only the signal:
				double DR = fSplineTimePDFDarkRate[pmtType]->Eval(residual);
				proba -= DR;
				// Then, scale signal so that signal integral / DR integral = signalOverNoise
				// To do so, we should scale signal so that integral = DR integral * signalOverNoise
				// And of course, to rescale, we should divide by the signal integral
				double Factor = fEventInfo.NoiseIntegral * fEventInfo.SignaloverNoise / fEventInfo.SignalIntegral;
				proba *= Factor;
				// And add again the DR
				proba += DR;
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3) std::cout << "proba after scaling=" << proba << std::endl;
#endif
			}
		}
		else if (killEdges)
		{
			if (residual >= upperLimit)
			{
				proba = fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3 && pmtType == 1)
					std::cout << "PMT type = " << pmtType << ", Upper limit = " << upperLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
			}
			else if (residual <= lowerLimit)
			{
				// continue;
				proba = fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3 && pmtType == 1)
					std::cout << "PMT type = " << pmtType << ", Lower limit = " << lowerLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
			}
		}
		else proba = 0;

		if (proba < 0)
		{
			std::cout << "Error in PDF" << std::endl;
			if (proba > -1e-1) proba = 0; // Since spline sometimes slightly goes below 0 due to interpolation
			else
			{
				std::cout << "Error in " << residual << "ns where proba = " << proba << std::endl;
				return 0;
			}
		}
		if (proba == 0) proba = 1e-20;
		NLL += -TMath::Log(proba);
	}
	
	//* Directionality, not used RN

// 	if (directionality != 0)
// 	{
// 		double NLLdir = 0;

// 		NLLdir = this->FindNLLDirectionality(vertexPosition, nhits, VERBOSE, lowerLimit, upperLimit);

// #ifdef VERBOSE_NLL
// 		if (VERBOSE >= 2)
// 			std::cout << "NLL = " << NLL << ", dir (L) = " << TMath::Exp(-NLLdir) << ", dir 2 = " << NLLdir << std::endl;
// #endif
// 		if (directionality == 1)
// 			NLL += NLLdir;
// 		else if (directionality == 2)
// 			NLL = NLLdir;
// 	}

	// timer.Stop();
	// std::cout << "NLL = " << NLL << std::endl;
	return NLL;
}

double Likelihoods::Vertex_Score(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double /*lowerLimit*/, double /*upperLimit*/, bool /*killEdges*/, bool /*scaleDR*/, int /*directionality*/)
{
	double NLL = 0;

	// std::cout << " Find NLL " << nhits << std::endl;
	// std::cout << " First Hit " << fHitCollection->At(0).PMT << " " << vertexPosition.size() << std::endl;
	for (int ihit = 0; ihit < nhits; ihit++)
	{
		// Hit lHit = fHitInfo[ihit];
		Hit lHit = lHitCol->At(ihit);

		// std::cout << " NLL Hit " << ihit << " " << lHit.PMT << std::endl;
		int iPMT = lHit.PMT;

		// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];

		int pmtType = Astro_GetPMTType(iPMT);
		// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

		// double tof = distance / fLightSpeed;
		// double residual = hitTime - tof - vertexPosition[3];
		double residual = ComputeResidualTime(vertexPosition, vertexPosition[3], lHitCol->At(ihit));
#ifdef KILLHALF
		if (pmtType == 1)
		{
			double t = fRand->Uniform(0, 1);
			if (t > 0.5)
				continue;
		}
#endif

		bool bCondition = (residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive) || (pmtType == 1 && !fLimit_mPMT);

		if (bCondition)
		{
			NLL++;
#ifdef VERBOSE_NLL
			if (VERBOSE >= 3)
			{
				std::cout << "Residual=" << residual << ", pmt type=" << pmtType << std::endl;
				std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << std::endl;
				std::cout << "residual=" << residual << std::endl;
			}
#endif
		}
	}

	return NLL != 0 ? -TMath::Log(NLL) : 1e15;
}

double Likelihoods::FindNLL(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality)
{
	double NLL = 0;

	for (int ihit = 0; ihit < nhits; ihit++)
	{
		Hit lHit = lHitCol->At(ihit);

		int iPMT = lHit.PMT;
		// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];

		int pmtType = Astro_GetPMTType(iPMT);
		// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

		// double tof = distance / fLightSpeed;
		// double residual = hitTime - tof - vertexPosition[3];
		double residual = ComputeResidualTime(vertexPosition, vertexPosition[3], lHitCol->At(ihit));
		double proba;
#ifdef KILLHALF
		if (pmtType == 1)
		{
			double t = fRand->Uniform(0, 1);
			if (t > 0.5)
				continue;
		}
#endif

		if (likelihood)
		{
			bool condition;

			condition = residual > lowerLimit && residual < upperLimit;
			if (condition)
			{
				proba = fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3)
				{
					std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << std::endl;
					std::cout << "Residual=" << residual << ", proba=" << proba << ", pmt type=" << pmtType << std::endl;
				}
#endif
				if (scaleDR)
				{
					std::cout << "scaleDR with Lower Limit: " << lowerLimit << " Upper Limit: " << upperLimit << std::endl;
					EventInfo fEventInfo = MakeEventInfo(lowerLimit, upperLimit, pmtType);

					// First, substract the DR from the PDF to keep only the signal:
					double DR = fSplineTimePDFDarkRate[pmtType]->Eval(residual);
					proba -= DR;
					// Then, scale signal so that signal integral / DR integral = signalOverNoise
					// To do so, we should scale signal so that integral = DR integral * signalOverNoise
					// And of course, to rescale, we should divide by the signal integral
					double Factor = fEventInfo.NoiseIntegral * fEventInfo.SignaloverNoise / fEventInfo.SignalIntegral;
					proba *= Factor;
					// And add again the DR
					proba += DR;
					if (VERBOSE >= 3)
						std::cout << "proba after scaling=" << proba << std::endl;
				}
			}
			else if (killEdges && residual >= upperLimit)
			{
				proba = fSplineTimePDFQueue[pmtType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3 && pmtType == 1)
					std::cout << "PMT type = " << pmtType << ", Upper limit = " << upperLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
			}
			else if (killEdges && residual <= lowerLimit)
			{
				proba = fSplineTimePDFQueue[pmtType]->Eval(lowerLimit);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3 && pmtType == 1)
					std::cout << "PMT type = " << pmtType << ", Lower limit = " << lowerLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
			}
			else
			{
				proba = 0;
			}
		}
		else
		{
			bool condition;
			if (pmtType == 0 || (pmtType == 1 && fLimit_mPMT)) condition = residual > fHitTimeLimitsNegative && residual < fHitTimeLimitsPositive;
			else condition = true;
			if (condition)
			{
				NLL++; //= fSplineTimePDFQueue[pmtType]->Eval(residual);
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3)
				{
					std::cout << "Residual=" << residual << ", pmt type=" << pmtType << std::endl;
					std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof=" << tof << std::endl;
					std::cout << "residual=" << residual << ", proba=" << proba << std::endl;
				}
#endif
			}
			else if (killEdges) continue;
			else proba = 0;
		}

		if (likelihood)
		{
			if (proba < 0)
			{
				std::cout << "Error in PDF" << std::endl;
				if (proba > -1e-1) proba = 0; // Since spline sometimes slightly goes below 0 due to interpolation
				else
				{
					std::cout << "Error in " << residual << "ns where proba = " << proba << std::endl;
					return 0;
				}
			}
			if (proba == 0) proba = 1e-20;
			NLL += -TMath::Log(proba);
		}
	}
	// std::cout << "FindNLL: Time it took for the loop = " << timer.RealTime() << std::endl;
	// cout<<"NLL="<<NLL<<endl;HitCollection
	if (!likelihood)
	{
		if (NLL != 0) NLL = -TMath::Log(NLL);
		else NLL = 1e15;
	}
	// if(VERBOSE>=2) cout<<"NLL="<<NLL<<endl;
	if (directionality != 0)
	{
		double NLLdir = 0;
		double NLLdir2 = 0;
		if (likelihood)
		{
			// NLLdir=findNLLDirectionalityBayes(vertexPosition, nhits, verbose,lowerLimit,upperLimit);
			NLLdir2 = FindNLLDirectionality(lHitCol, vertexPosition, nhits, verbose, lowerLimit, upperLimit);

			// double test = FindNLLDirectionality(lHitCol, fHitCollection, vertexPosition, nhits, verbose, lowerLimit, upperLimit);
			// std::cout << "test = " << test << std::endl;
		}
		else
		{
			NLLdir2 = FindNLLDirectionality(lHitCol, vertexPosition, nhits, verbose, fHitTimeLimitsNegative, fHitTimeLimitsPositive);
		}
		if (VERBOSE >= 2) std::cout << "NLL = " << NLL << ", dir (L) = " << TMath::Exp(-NLLdir) << ", dir 2 = " << NLLdir2 << std::endl;
		if (directionality == 1) NLL += NLLdir2;
		else if (directionality == 2) NLL = NLLdir2;
	}

	return NLL;
}

double Likelihoods::ComputeDirNLL_NoPDF(const HitCollection<Hit>* lHitCol, const std::vector<double>& vertexPosition, const std::vector<double>& vertexDirection, int nhits) {
	double DNLL = 0.;
    for (int ihit = 0; ihit < nhits; ihit++) 
	{
        Hit lHit = lHitCol->At(ihit);
        int iPMT = lHit.PMT;

        PMTInfo lPMTInfo = (*fPMTList)[iPMT];

        // Calculate vector from PMT to vertex
        double dirVector[3];
        dirVector[0] = -(vertexPosition[0] - lPMTInfo.Position[0]);
        dirVector[1] = -(vertexPosition[1] - lPMTInfo.Position[1]);
        dirVector[2] = -(vertexPosition[2] - lPMTInfo.Position[2]);
		double vertexdirnorm= pow(vertexDirection[0]*vertexDirection[0]+vertexDirection[1]*vertexDirection[1]+vertexDirection[2]*vertexDirection[2],0.5);
		double vertexDirNormalized[3];

		for (int j = 0; j<3;j++) vertexDirNormalized[j]=vertexDirection[j]/vertexdirnorm;

        // Normalize the vector
        double length = sqrt(dirVector[0] * dirVector[0] + dirVector[1] * dirVector[1] + dirVector[2] * dirVector[2]);
        dirVector[0] /= length;
        dirVector[1] /= length;
        dirVector[2] /= length;
		
        // Calculate RelativeAngle (cosine of angle between vectors)
        double RelativeAngle = dirVector[0] * vertexDirNormalized[0] + dirVector[1] * vertexDirNormalized[1] + dirVector[2] * vertexDirNormalized[2];
		double theta = acos(RelativeAngle)*180./TMath::Pi();
		//std::cout << "theta: " << theta << std::endl;
        if (theta<46 && theta>38)
		{
			std::cout << "Got hit! " << theta << std::endl;
        	DNLL += -1;
		} 

		
		if (VERBOSE >= 2) std::cout << "Direction NLL = " << DNLL << std::endl;
	}
	
    return DNLL;

}

//Uses the PDF to find the DirNLL, used mainly for the MIGRAD optimization
double Likelihoods::Dir_NLL(const HitCollection<Hit>* lHitCol, const std::vector<double>& vertexPosition, double theta_track, double phi_track, int nhits) {
    double DNLL = 0;

    // Ensure theta and phi are within valid ranges
    // theta in [0, pi], phi in [-pi, pi]
    if (theta_track < 0 || theta_track > TMath::Pi()) 
	{
        std::cerr << "Error: theta_track is out of range [0, π]." << std::endl;
        return std::numeric_limits<double>::infinity();
    }
    if (phi_track < -TMath::Pi() || phi_track > TMath::Pi()) 
	{
        std::cerr << "Error: phi_track is out of range [-π, π]." << std::endl;
        return std::numeric_limits<double>::infinity();
    }

    // Precompute sine and cosine of theta_track and phi_track
    // double sin_theta_track = sin(theta_track);
    // double cos_theta_track = cos(theta_track);
    // double sin_phi_track = sin(phi_track);
    // double cos_phi_track = cos(phi_track);
	/* 1
	double a = 0.556805;
	double b = 1.06697;
	double c = -0.902084;
	double d = -0.0707118;*/

	// 2
	// double a = -0.000440808;
	// double b = 1.71313;
	// double c = -2.01675;
	// double d = 1.81755;
    
	// Direction vector components from theta and phi
    // double vertexDirNormalized[3];
    // vertexDirNormalized[0] = sin_theta_track * cos_phi_track;
    // vertexDirNormalized[1] = sin_theta_track * sin_phi_track;
    // vertexDirNormalized[2] = cos_theta_track;

    for (int ihit = 0; ihit < nhits; ihit++) {
		int trueNHits = lHitCol->Size();
		if(ihit >= trueNHits)
		{
			std::cout << "Error: ihit is out of range in direction NLL computation" << std::endl;
			return std::numeric_limits<double>::infinity();
		}
		if(vertexPosition.size() != 4)
		{
			std::cout << "Error: vertexPosition is not of size 4 in direction NLL computation" << std::endl;
			return std::numeric_limits<double>::infinity();
		}

		Hit lHit = lHitCol->At(ihit);
		int iPMT = lHit.PMT;

		// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
		PMTInfo lPMTInfo = (*fPMTList)[iPMT];
		
		int pmtType = Astro_GetPMTType(iPMT);
		// double distance = Astro_GetDistance(lPMTInfo.Position,vertexPosition);
		
		// double tof = distance / fLightSpeed;
		// double residual = hitTime - tof - vertexPosition[3];

		double residual = ComputeResidualTime(vertexPosition, vertexPosition[3], lHitCol->At(ihit));
	
		bool DirCondition = (residual >= -5 && residual <= 15) || (pmtType == 1 && 	!fLimit_mPMT);


		if(DirCondition || DirTakeAll)
		{		
			int iPMT = lHit.PMT;
			PMTInfo lPMTInfo = (*fPMTList)[iPMT];

			// Ensure lPMTInfo.Position and lPMTInfo.Orientation have the correct size
			if (sizeof(lPMTInfo.Position) / sizeof(lPMTInfo.Position[0]) < 3 || 
				sizeof(lPMTInfo.Orientation) / sizeof(lPMTInfo.Orientation[0]) < 3) {
				std::cerr << "Error: lPMTInfo.Position or lPMTInfo.Orientation is not of size 3." << std::endl;
				continue;
			}

			double dirVector[3];
			dirVector[0] = lPMTInfo.Position[0] - vertexPosition[0];
			dirVector[1] = lPMTInfo.Position[1] - vertexPosition[1];
			dirVector[2] = lPMTInfo.Position[2] - vertexPosition[2];

			double dirNorm = sqrt(dirVector[0]*dirVector[0] + dirVector[1]*dirVector[1] + dirVector[2]*dirVector[2]);
			if (dirNorm == 0)
			{
				std::cerr << "Error: dirNorm is zero, cannot normalize direction vector." << std::endl;
				continue;
			}
			dirVector[0] /= dirNorm;
			dirVector[1] /= dirNorm;
			dirVector[2] /= dirNorm;

			// Convert dirVector to spherical coordinates (theta_pmt, phi_pmt)fDirectionPDF
			double theta_pmt = acos(dirVector[2]); // theta in [0, π]
			double phi_pmt = atan2(dirVector[1], dirVector[0]); // phi in [-π, π]

			// Calculate the angular difference between the two directions
			// double delta_theta = theta_pmt - theta_track; //* could be useful
			double delta_phi = phi_pmt - phi_track;

			// Ensure delta_phi is within [-π, π]
			if (delta_phi > TMath::Pi()) delta_phi -= 2 * TMath::Pi();
			else if (delta_phi < -TMath::Pi()) delta_phi += 2 * TMath::Pi();

			// Calculate the angle between the two directions using the spherical law of cosines
			double cos_relative_angle = sin(theta_track) * sin(theta_pmt) * cos(delta_phi) +
										cos(theta_track) * cos(theta_pmt);

			//cos_relative_angle = std::max(-1.0, std::min(1.0, cos_relative_angle));

			double relative_angle = acos(cos_relative_angle); // Angle in radians

			// Convert angle to degrees if your PDF is defined in degrees
			double theta = relative_angle * 180.0 / TMath::Pi();

			// Evaluate your PDF at theta
			// std::cout << "before pdf : " << ihit << std::endl;
			double proba = fDirectionPDF->Eval(cos_relative_angle)/*sin(theta)*/;
			// std::cout << "after pdf : " << ihit << std::endl;

			// Correct for PMT orientation if necessary
			double PMTOrientation[3];
			PMTOrientation[0] = lPMTInfo.Orientation[0];
			PMTOrientation[1] = lPMTInfo.Orientation[1];
			PMTOrientation[2] = lPMTInfo.Orientation[2];

			// Normalize PMT orientation vector
			double pmtOrientNorm = sqrt(PMTOrientation[0]*PMTOrientation[0] +
										PMTOrientation[1]*PMTOrientation[1] +
										PMTOrientation[2]*PMTOrientation[2]);

			PMTOrientation[0] /= pmtOrientNorm;
			PMTOrientation[1] /= pmtOrientNorm;
			PMTOrientation[2] /= pmtOrientNorm;

			double cosPMThitAngle = dirVector[0]*PMTOrientation[0] +
											dirVector[1]*PMTOrientation[1] +
											dirVector[2]*PMTOrientation[2];
			// double PMThitAngle = acos(-cosPMThitAngle);
			cosPMThitAngle *= -1;
			// double cos2PMThitAngle = fabs(-2 * cosPMThitAngle * cosPMThitAngle + 1);
			// double Ftheta = a + b*cosPMThitAngle + c*(pow(cosPMThitAngle,2)) + d*(pow(cosPMThitAngle,3));
			
			//proba *= hitCharge;

			// Avoid log(0) by ensuring proba is positive
			if (proba <= 0) 
			{
				proba = 1e-20;
				if (VERBOSE >= 1) std::cout << "Warning: Probability is zero or negative at hit #" << ihit << std::endl;
			}
			//std::cout << "Proba Corrected by " << PMThitAngle * 180 / M_PI << " degrees with a sincos " << sin(PMThitAngle) << ", " << cosPMThitAngle << std::endl;
			// Accumulate Negative Log-Likelihood
			DNLL += -TMath::Log(proba);//*(cosPMThitAngle/Ftheta)*(sin(PMThitAngle));
		
			if (VERBOSE >= 3) std::cout << "Hit #" << ihit << ", Theta = " << theta << " degrees, Probability = " << proba << std::endl;
    	}

		if (VERBOSE >= 2) std::cout << "Total Direction NLL = " << DNLL << std::endl;
	}
    return DNLL;
}


double Likelihoods::FindNLLDirectionality(const HitCollection<Hit>* lHitCol, std::vector<double> vVtxPos, int nhits, int verbose, double /*lowerLimit*/, double /*upperLimit*/)
{

	double NLL = 0;
	// std::vector<double> vDirection(2, 0.); // Return phi and theta.

	// // Interpolate isn't compatible with multi-thread, need to use mutex which lead to long deadtime.
	// // Copy TGraph2D
	// mtx.lock();
	// static thread_local TGraph2D gPMTDirectionality_2D_local_0 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][0]);
	// static thread_local TGraph2D gPMTDirectionality_2D_local_1 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][1]);
	// static thread_local TGraph2D gPMTDirectionality_2D_local_2 = TGraph2D(*gPMTDirectionality_2D[MiniPMT][2]);
	// mtx.unlock();

	// TGraph2D *tDirectionality[3];
	// tDirectionality[0] = &gPMTDirectionality_2D_local_0;
	// tDirectionality[1] = &gPMTDirectionality_2D_local_1;
	// tDirectionality[2] = &gPMTDirectionality_2D_local_2;

	// for (int ihit = 0; ihit < nhits; ihit++)
	// {
	// 	// Hit lHit = fHitInfo[ihit];
	// 	Hit lHit = LeafInputs::fHitCollection->At(ihit);

	// 	int iPMT = lHit.PMT;

	// 	PMTInfo lPMTInfo = (*fPMTList)[iPMT];
	// 	int pmtType = Astro_GetPMTType(iPMT);

	// 	if (pmtType == 0)
	// 		continue;
	// 	bool condition = true;

	// 	if (condition)
	// 	{
	// 		double vPMTVtx[4];
	// 		this->VectorVertexPMT(vVtxPos, iPMT, vPMTVtx);

	// 		double dPhi = vPMTVtx[0];
	// 		double dTheta = vPMTVtx[1];
	// 		double dDist = vPMTVtx[2];

	// 		int pmtGroup = lPMTInfo.mPMT_Group;

	// 		double proba = 0;
	// 		proba = tDirectionality[pmtGroup]->Interpolate(dPhi, dTheta);

	// 		double dDistCorr = fDistResponsePMT[pmtType]->Eval(dDist);

	// 		proba *= dDistCorr;

	// 		if (proba == 0)
	// 		{
	// 			proba = 1e-20;
	// 			if (verbose)
	// 				std::cout << "We are at proba = 0, theta = " << dTheta << ", proba used = " << proba << std::endl;
	// 		}
	// 		NLL += -TMath::Log(proba);

	// 		if (VERBOSE >= 3)
	// 		{
	// 			// std::cout<<"PMT type="<<pmtType<< ", hit#"<<ihit<<", theta =" << vDirection[1] << ", proba="<<proba<<std::endl;
	// 			std::cout << "PMT type=" << pmtType << ", hit#" << ihit << ", theta =" << dTheta << ", proba=" << proba << std::endl;
	// 		}
	// 	}
	// }
	// if (VERBOSE >= 2)
	// 	std::cout << "NLL directionnel=" << NLL << std::endl;
	return NLL;
}

double Likelihoods::GoodnessOfFit(const HitCollection<Hit>* lHitCol, std::vector<double> vertexPosition, int nhits, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality)
{
	if (!lHitCol || nhits == 0) return -1;

    // double NLL_theory = Vertex_Time_NLL(lHitCol, vertexPosition, nhits, lowerLimit, upperLimit, killEdges, scaleDR, directionality);

	double NLL_theory = vertexPosition[4];

    std::vector<double> residuals;
	for(int ihit = 0; ihit < nhits; ihit++) residuals.push_back(ComputeResidualTime(vertexPosition, vertexPosition[3], lHitCol->At(ihit)));

    double bandwidth = EstimateBandeWidth(residuals);

    //* Compute KDE and total log likelihood of data
    double log_p_data = 0.0;
    for (const auto& t_i : residuals) 
	{
        double p = KDE_Estimate(residuals, t_i, bandwidth);
        log_p_data += std::log(p + 1e-12);  // avoid log(0)
    }

    //* Likelihood ratio
    double NLL_data = -log_p_data;
    double NLLR = NLL_data - NLL_theory;

    return NLLR;
}