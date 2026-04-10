/*****************************************************************************************************/
/**	LEAFSrcInputs.cc																				**/
/**	Date: March 27th 2026																			**/
/**	Desc: Implementation of the likelihood functions												**/
/*****************************************************************************************************/

#ifndef LEAF_HPP
#include "LEAF.hpp"
#endif

//Function to optimize, for MIGRAD
void LEAFLikelihoods::MinuitLikelihood(int & /*nDim*/, double * /*gout*/, double &NLL, double par[], int /*flg*/)
{
	ROOT::Math::XYZTVector vertex(par[0], par[1], par[2], par[3]); // In centimeters
	ROOT::Math::XYZVector direction(par[10], par[11], par[12]);

	// double dirThreshold = par[13];
	double lowerLimit = par[6];
	double upperLimit = par[7];
	double directionality = par[9];
	// std::vector<double> vertexDirection(3,0.);
	double timeNLL = LEAF::GetME()->Vertex_Time_NLL(vertex, lowerLimit, upperLimit, true, false, directionality);
	// double angleNLL = LEAF::GetME()->AngleNLL(vertex, fTrueDir);
	NLL = timeNLL;
}

//Function to optimize, for MIGRAD
void LEAFLikelihoods::MinuitDirNLL(int& /*nDim*/, double* /*gout*/, double& DNLL, double par[], int /*flg*/) 
{
	// Extract theta and phi from parameters
	double theta = par[0];
	double phi = par[1];

	// Ensure the vertex is properly initialized
	ROOT::Math::XYZTVector vertex(par[3], par[4], par[5], par[6]); 

	// Ensure theta and phi are within valid ranges
	if (theta < 0 || theta > TMath::Pi()) {
		DNLL = 1e10; // Assign a large NLL to invalid directions
		return;
	}
	if (phi < -TMath::Pi() || phi > TMath::Pi()) {
		DNLL = 1e10;
		return;
	}

	// Calculate NLL using theta and phi
	DNLL = LEAF::GetME()->Dir_NLL(vertex, theta, phi);
}

void LEAFLikelihoods::MinuitJointNLL(int& nDim, double * gout, double & NLL, double par[], int flg)
{
	ROOT::Math::XYZTVector vertex(par[0], par[1], par[2], par[3]);
	double theta = par[4];
	double phi = par[5];
	double lowerLimit = par[8];
	double upperLimit = par[9];
	// double directionality = par[9];

	// Ensure theta and phi are within valid ranges
	if (theta < 0 || theta > TMath::Pi()) {
		NLL = 1e10; // Assign a large NLL to invalid directions
		return;
	}
	if (phi < -TMath::Pi() || phi > TMath::Pi()) {
		NLL = 1e10;
		return;
	}

	double VtxNLL = LEAF::GetME()->Vertex_Time_NLL(vertex, lowerLimit, upperLimit, true, false, 0);
	double DirNLL = LEAF::GetME()->Dir_NLL(vertex, theta, phi);

	NLL = VtxNLL * DirNLL;
}

double LEAF::Vertex_Time_NLL(ROOT::Math::XYZTVector vertex, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality)
{
	double NLL = 0;

	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());

	this->GenerateEventInfo(lowerLimit, upperLimit);

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);

		// Loop over ID PMT hits
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);

			int iPMT = lHit->GetPMTNumber();
			// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

			// double tof = distance / fLightSpeed;
			// double residual = hitTime - tof - vertexPosition[3];
			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, lHit->GetTime(), vertex.T());

			// std::vector<double> toHit = VectorToHitNorm(vertexPosition, lHit);
			// if(dot(fTrueDir, toHit) < cos(90 * TMath::DegToRad()))
			// if(dot(fTrueDir, toHit) < 0)
			// {
			// 	continue; //! Skip this hit if it's not in the same hemisphere
			// }

			double proba = 0;

			bool condition = residual > lowerLimit && residual < upperLimit;

			if (condition)
			{
				proba = fSplineTimePDFQueue[iPMTType]->Eval(residual);

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

					// First, substract the DR from the PDF to keep only the signal:
					double DR = fSplineTimePDFDarkRate[iPMTType]->Eval(residual);
					proba -= DR;
					// Then, scale signal so that signal integral / DR integral = signalOverNoise
					// To do so, we should scale signal so that integral = DR integral * signalOverNoise
					// And of course, to rescale, we should divide by the signal integral
					double Factor = fEventInfo[iPMTType].NoiseIntegral * fEventInfo[iPMTType].SignalOverNoise / fEventInfo[iPMTType].SignalIntegral;
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
					proba = fSplineTimePDFQueue[iPMTType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
					if (VERBOSE >= 3 && pmtType == PMTType::kmPMT)
						std::cout << "PMT type = " << pmtType << ", Upper limit = " << upperLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
				}
				else if (residual <= lowerLimit)
				{
					// continue;
					proba = fSplineTimePDFQueue[iPMTType]->Eval(lowerLimit);
#ifdef VERBOSE_NLL
					if (VERBOSE >= 3 && pmtType == PMTType::kmPMT)
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

	// 		NLLdir = this->FindNLLDirectionality(vertexPosition, VERBOSE, lowerLimit, upperLimit);

	// #ifdef VERBOSE_NLL
	// 		if (VERBOSE >= 2)
	// 			std::cout << "NLL = " << NLL << ", dir (L) = " << TMath::Exp(-NLLdir) << ", dir 2 = " << NLLdir << std::endl;
	// #endif
	// 		if (directionality == 1)
	// 			NLL += NLLdir;
	// 		else if (directionality == 2)
	// 			NLL = NLLdir;
	// 	}

	}
	//timer.Stop();
	//std::cout << "NLL = " << NLL << std::endl;
	return NLL;
}

double LEAF::Vertex_Score(ROOT::Math::XYZTVector vertex, double /*lowerLimit*/, double /*upperLimit*/, bool /*killEdges*/, bool /*scaleDR*/, int /*directionality*/)
{
	double NLL = 0;
	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());

	// std::cout << " First Hit " << fHitCollection->At(0).PMT << " " << vertexPosition.size() << std::endl;
	
	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			// Hit lHit = fHitInfo[ihit];
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);

			int iPMT = lHit->GetPMTNumber();
			// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

			// double tof = distance / fLightSpeed;
			// double residual = hitTime - tof - vertexPosition[3];
			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, lHit->GetTime(), vertex.T());
			
			/*
			if (ihit < 5) {
				double distance = LEAFUtilities::GetDistance(vertexPos,lPMTInfo->GetPositionInCm());
				std::cout << "PMT: " << iPMT << " \t X: " << lPMTInfo->GetPositionInCm().X() << " Y: " << lPMTInfo->GetPositionInCm().Y() << " Z: " << lPMTInfo->GetPositionInCm().Z() << " T " << lHit->GetTime() << std::endl;
				std::cout << "Vertex: \t X:" << vertexPos.X() << " Y: " << vertexPos.Y() << " Z: " << vertexPos.Z() << " T " << vertex.T() << std::endl;
				ROOT::Math::XYZVector distVector = (vertexPos - lPMTInfo->GetPositionInCm());
				std::cout << "Diff: \t X:" << distVector.X() << " Y: " << distVector.Y() << " Z: " << distVector.Z() << std::endl;
				std::cout << " Distance " << distance << " -> " << residual << std::endl;
			}
			*/
#ifdef KILLHALF
			if (pmtType == PMTType::kmPMT) {
				double t = fRand->Uniform(0, 1);
				if (t > 0.5)
					continue;
			}
#endif

			bool bCondition = (residual > LEAFConfig::fHitTimeLimitsNegative && residual < LEAFConfig::fHitTimeLimitsPositive) || (pmtType == PMTType::kmPMT && !LEAFConfig::fLimit_mPMT);

			if (bCondition) {
				NLL++;
#ifdef VERBOSE_NLL
				if (VERBOSE >= 3) {
					std::cout << "Residual=" << residual << ", pmt type=" << pmtType << std::endl;
					std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << std::endl;
					std::cout << "residual=" << residual << std::endl;
				}
#endif
			}
		}
	}

	return NLL != 0 ? -TMath::Log(NLL) : 1e15;
}

double LEAF::FindNLL(ROOT::Math::XYZTVector vertex, bool likelihood, int verbose, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality)
{
	double NLL = 0;
	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());

	this->GenerateEventInfo(lowerLimit, upperLimit);

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
	
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
			// double distance = Astro_GetDistance(lPMTInfo.Position, vertexPosition);

			// double tof = distance / fLightSpeed;
			// double residual = hitTime - tof - vertexPosition[3];
			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, lHit->GetTime(), vertex.T());
			double proba;
	#ifdef KILLHALF
			if (pmtType == PMTType::kmPMT)
			{
				double t = fRand->Uniform(0, 1);
				if (t > 0.5)
					continue;
			}
	#endif

			if (likelihood) {
				bool condition;

				condition = residual > lowerLimit && residual < upperLimit;
				if (condition) {
					proba = fSplineTimePDFQueue[iPMTType]->Eval(residual);
#ifdef VERBOSE_NLL
					if (VERBOSE >= 3) {
						std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << std::endl;
						std::cout << "Residual=" << residual << ", proba=" << proba << ", pmt type=" << pmtType << std::endl;
					}
#endif
					if (scaleDR) {
						std::cout << "scaleDR with Lower Limit: " << lowerLimit << " Upper Limit: " << upperLimit << std::endl;

						// First, substract the DR from the PDF to keep only the signal:
						double DR = fSplineTimePDFDarkRate[iPMTType]->Eval(residual);
						proba -= DR;
						// Then, scale signal so that signal integral / DR integral = signalOverNoise
						// To do so, we should scale signal so that integral = DR integral * signalOverNoise
						// And of course, to rescale, we should divide by the signal integral
						double Factor = fEventInfo[iPMTType].NoiseIntegral * fEventInfo[iPMTType].SignalOverNoise / fEventInfo[iPMTType].SignalIntegral;
						proba *= Factor;
						// And add again the DR
						proba += DR;
						if (VERBOSE >= 3)
							std::cout << "proba after scaling=" << proba << std::endl;
					}
				}
				else if (killEdges && residual >= upperLimit) {
					proba = fSplineTimePDFQueue[iPMTType]->Eval(upperLimit);
#ifdef VERBOSE_NLL
					if (VERBOSE >= 3 && pmtType == PMTType::kmPMT)
						std::cout << "PMT type = " << pmtType << ", Upper limit = " << upperLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
#endif
				}
				else if (killEdges && residual <= lowerLimit) {
					proba = fSplineTimePDFQueue[iPMTType]->Eval(lowerLimit);
	#ifdef VERBOSE_NLL
					if (VERBOSE >= 3 && pmtType == PMTType::kmPMT)
						std::cout << "PMT type = " << pmtType << ", Lower limit = " << lowerLimit << ", residual = " << residual << ", proba = " << proba << std::endl;
	#endif
				}
				else {
					proba = 0;
				}
			}
			else {
				bool condition;
				if (pmtType == PMTType::kID || (pmtType == PMTType::kmPMT && LEAFConfig::fLimit_mPMT)) {
					condition = (residual > LEAFConfig::fHitTimeLimitsNegative && residual < LEAFConfig::fHitTimeLimitsPositive);
				}
				else {
					condition = true;
				}
				if (condition) {
					NLL++; //= fSplineTimePDFQueue[pmtType]->Eval(residual);
	#ifdef VERBOSE_NLL
					if (VERBOSE >= 3) {
						std::cout << "Residual=" << residual << ", pmt type=" << pmtType << std::endl;
						std::cout << "hit#" << ihit << ", hit time =" << hitTime << ", vertex time = " << vertexPosition[3] << ", distance PMT vs vertex = " << distance << ", tof=" << tof << std::endl;
						std::cout << "residual=" << residual << ", proba=" << proba << std::endl;
					}
	#endif
				}
				else if (killEdges) {
					continue;
				}
				else {
					proba = 0;
				}
			}

			if (likelihood) {
				if (proba < 0) {
					std::cout << "Error in PDF" << std::endl;
					if (proba > -1e-1) proba = 0; // Since spline sometimes slightly goes below 0 due to interpolation
					else {
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
		if (!likelihood) {
			if (NLL != 0) {
				NLL = -TMath::Log(NLL);
			}
			else {
				NLL = 1e15;
			}
		}
		// if(VERBOSE>=2) cout<<"NLL="<<NLL<<endl;
		if (directionality != 0) {
			double NLLdir = 0;
			double NLLdir2 = 0;
			if (likelihood) {
				// NLLdir=findNLLDirectionalityBayes(vertex, verbose,lowerLimit,upperLimit);
				NLLdir2 = this->FindNLLDirectionality(vertex, verbose, lowerLimit, upperLimit);

				// double test = FindNLLDirectionality(lHitCol, fHitCollection, vertexPosition, verbose, lowerLimit, upperLimit);
				// std::cout << "test = " << test << std::endl;
			}
			else {
				NLLdir2 = this->FindNLLDirectionality(vertex, verbose, LEAFConfig::fHitTimeLimitsNegative, LEAFConfig::fHitTimeLimitsPositive);
			}
			if (VERBOSE >= 2) {
				std::cout << "NLL = " << NLL << ", dir (L) = " << TMath::Exp(-NLLdir) << ", dir 2 = " << NLLdir2 << std::endl;
			}

			if (directionality == 1) {
				NLL += NLLdir2;
			}
			else if (directionality == 2) {
				NLL = NLLdir2;
			}
		}
	}

	return NLL;
}

double LEAF::ComputeDirNLL_NoPDF(const ROOT::Math::XYZTVector& vertex, const ROOT::Math::XYZVector& direction) {
	double DNLL = 0.;

	ROOT::Math::XYZVector vertexPos(vertex.X(),vertex.Y(),vertex.Z());

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			// Calculate vector from PMT to vertex
			ROOT::Math::XYZVector dirVectorNormalized = (vertexPos - lPMTInfo->GetPositionInCm()).Unit();
			ROOT::Math::XYZVector directionNormalized = direction.Unit();
			
			// Calculate RelativeAngle (cosine of angle between vectors)
			double relativeAngle = dirVectorNormalized.Dot(directionNormalized);
			double theta = acos(relativeAngle)*180./TMath::Pi();
			//std::cout << "theta: " << theta << std::endl;
			if (theta<46 && theta>38) {
				std::cout << "Got hit! " << theta << std::endl;
				DNLL += -1;
			} 
						
			if (VERBOSE >= 2) std::cout << "Direction NLL = " << DNLL << std::endl;
		}
	}
	
    return DNLL;

}

//Uses the PDF to find the DirNLL, used mainly for the MIGRAD optimization
double LEAF::Dir_NLL(const ROOT::Math::XYZTVector& vertex, double theta_track, double phi_track) 
{
    double DNLL = 0;
	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());

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

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			// double hitTime = (fTimeCorrection + lHit.T) / TimeDelta::ns;
			
			// double distance = Astro_GetDistance(lPMTInfo.Position,vertexPosition);
			
			// double tof = distance / fLightSpeed;
			// double residual = hitTime - tof - vertexPosition[3];

			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, lHit->GetTime(), vertex.T());

			bool DirCondition = (residual >= LEAFConfig::fDirectionPDF_minResidual && residual <= LEAFConfig::fDirectionPDF_maxResidual) || (pmtType == PMTType::kmPMT && !LEAFConfig::fLimit_mPMT);


			if(DirCondition || LEAFConfig::fDirTakeAll) {		

				ROOT::Math::XYZVector dirVector = (lPMTInfo->GetPositionInCm() - vertexPos).Unit();
				
				// Convert dirVector to spherical coordinates (theta_pmt, phi_pmt)fDirectionPDF
				double theta_pmt = acos(dirVector.Z()); // theta in [0, π]
				double phi_pmt = atan2(dirVector.Y(), dirVector.X()); // phi in [-π, π]

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
				double proba = fDirectionPDF[iPMTType]->Eval(theta); //* Nicolas

				// Correct for PMT orientation if necessary
				ROOT::Math::XYZVector pmtOrientation = lPMTInfo->GetOrientation().Unit();

				double cosPMThitAngle = dirVector.Dot(pmtOrientation);
				// double PMThitAngle = acos(-cosPMThitAngle);
				cosPMThitAngle *= -1;
				// double cos2PMThitAngle = fabs(-2 * cosPMThitAngle * cosPMThitAngle + 1);
				// double Ftheta = a + b*cosPMThitAngle + c*(pow(cosPMThitAngle,2)) + d*(pow(cosPMThitAngle,3));
				
				//proba *= hitCharge;

				// Avoid log(0) by ensuring proba is positive
				if (proba <= 0) {
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
	}
    return DNLL;
}

double LEAF::AngleNLL(ROOT::Math::XYZTVector vertex, ROOT::Math::XYZVector direction)
{
	double NLL = 0;
	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());
	
	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			//* compute angle
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);
			
			ROOT::Math::XYZVector toPMT = (lPMTInfo->GetPositionInCm() - vertexPos).Unit();
			double dotP = toPMT.Dot(direction.Unit());
			
			dotP = std::max(-1.0, std::min(1.0, dotP)); // Ensure dot product is within valid range for acos
			double angle = std::acos(dotP) * 180 / TMath::Pi();
			NLL += -TMath::Log(std::max(1e-20, fDirectionPDF[iPMTType]->Eval(angle))); //? Nicolas
		}
	}

	return NLL;
}


double LEAF::FindNLLDirectionality(ROOT::Math::XYZTVector vVtxPos, int verbose, double /*lowerLimit*/, double /*upperLimit*/) {

	double NLL = 0;
	return NLL;
}


double LEAF::GoodnessOfFit(ROOT::Math::XYZTVector vertex, double NLL_theory, double lowerLimit, double upperLimit, bool killEdges, bool scaleDR, int directionality) {
    // double NLL_theory = Vertex_Time_NLL(lHitCol, vertexPosition, lowerLimit, upperLimit, killEdges, scaleDR, directionality);

	ROOT::Math::XYZVector vertexPos(vertex.X(), vertex.Y(), vertex.Z());
	std::vector<double> residuals;

	for (auto pmtType : LEAFConfig::fActivePMTTypes) {
		const int iPMTType = static_cast<int>(pmtType);
		for (unsigned int ihit =0; ihit < fHitsCollection[iPMTType]->size(); ihit++) {
			const HKHit* lHit = fHitsCollection[iPMTType]->at(ihit);
			int iPMT = lHit->GetPMTNumber();
			const HKGeometryPMT* lPMTInfo = fGeoPMTs[iPMTType]->at(iPMT);

			double residual = LEAFUtilities::GetResidual(lPMTInfo->GetPositionInCm(), vertexPos, lHit->GetTime(), vertex.T());
			residuals.push_back(residual);
		}
	}

    double bandwidth = LEAFUtilities::EstimateBandeWidth(residuals);

    //* Compute KDE and total log likelihood of data
    double log_p_data = 0.0;
    for (const auto& t_i : residuals) {
        double p = LEAFUtilities::KDE_Estimate(residuals, t_i, bandwidth);
        log_p_data += std::log(p + 1e-12);  // avoid log(0)
    }

    //* Likelihood ratio
    double NLL_data = -log_p_data;
    double NLLR = NLL_data - NLL_theory;

    return NLLR;
}