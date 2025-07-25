#include "LeafUtility.hh"

bool ContainsTrueVtx(std::vector<VtxCandidate>* candidates, std::vector<double> fTrueVtxPos)
{
	for (unsigned int i = 0; i < candidates->size(); i++)
	{
		std::vector<double> vtx = {(*candidates)[i].X, (*candidates)[i].Y, (*candidates)[i].Z, (*candidates)[i].T};
		double distance = Distance3D(vtx, fTrueVtxPos);
		if (distance <= fSearchVtxStep) //? should we take step / 2 ?
		{
			return true;
		}
	}
	return false;
}

bool ContainsTrueDir(std::vector<DirectionCandidate>* candidates, std::vector<double> fTrueDir)
{
	for (unsigned int i = 0; i < candidates->size(); i++)
	{
		std::vector<double> truePolarDir = CartesianToPolarNorm(fTrueDir);
		if(abs((*candidates)[i].theta - truePolarDir[0]) < theta_step && abs((*candidates)[i].phi - truePolarDir[1]) < phi_step)
		{
			return true;
		}
	}
	return false;
}


//* to know if the predicted vertex and the true vertex belong to the same candidate
//? doesn't work because the radius overlap between candidates
bool CanidatesMatch(std::vector<std::vector<double>>* candidates, std::vector<double> point, std::vector<double> fTrueVtxPos)
{
	for(unsigned int i = 0; i < candidates->size(); i++)
	{
		if (Distance3D((*candidates)[i], fTrueVtxPos) <= fSearchVtxStep && Distance3D((*candidates)[i], point) <= fSearchVtxStep) return true;
	}
	return false;
}

double Distance3D(std::vector<double> point1, std::vector<double> point2)
{
	return sqrt(pow(point1[0] - point2[0], 2) + 
				pow(point1[1] - point2[1], 2) + 
				pow(point1[2] - point2[2], 2));
}

void Normalize(double a[3])
{
	double l;
	l=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	for(int j=0; j<3; j++){
		a[j]/=l;
	}
}

std::vector<double>& Normalize(std::vector<double>& vector)
{
	if (vector.size() < 3) throw std::invalid_argument("Vector must be of size at least 3 to be normalized");
	double length = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
	if (length == 0) throw std::invalid_argument("Cannot normalize a zero-length vector.");
	for (int i = 0; i < 3; ++i) vector[i] /= length;
	return vector;
}

double calculateDistance(const std::vector<double>& A, const std::vector<double>& B) 
{
    return sqrt(pow(A[0] - B[0], 2) + pow(A[1] - B[1], 2) + pow(A[2] - B[2], 2));
}

double dot(const std::vector<double>& A, const std::vector<double>& B)
{
	if (A.size() != B.size()) throw std::invalid_argument("Vectors must be of the same size for dot product.");
	double result = 0;
	for (size_t i = 0; i < A.size(); ++i) result += A[i] * B[i];
	return result;
}

std::vector<double> PolarToCartesianNorm(std::vector<double> polarVector)
{
	double theta = polarVector[0];
	double phi = polarVector[1];
	std::vector<double> cartVector = std::vector<double>(3);
	cartVector[0] = sin(theta) * cos(phi);
	cartVector[1] = sin(theta) * sin(phi);
	cartVector[2] = cos(theta);
	Normalize(cartVector);
	return cartVector;
}

std::vector<double> CartesianToPolarNorm(const std::vector<double>& cartVector)
{
    std::vector<double> polarVector(2);
    double x = cartVector[0];
    double y = cartVector[1];
    double z = cartVector[2];

    polarVector[0] = acos(z);              // theta
    polarVector[1] = atan2(y, x);          // phi

    return polarVector;
}

double EstimateBandeWidth(const std::vector<double>& residuals) 
{
	double sum = std::accumulate(residuals.begin(), residuals.end(), 0.0);
	double mean = sum / residuals.size();

	double sq_sum = std::inner_product(residuals.begin(), residuals.end(), residuals.begin(), 0.0);
	double variance = sq_sum / residuals.size() - mean * mean;
	double stddev = std::sqrt(variance);
	return 1.06 * stddev * std::pow(residuals.size(), -1.0 / 5.0);
}

double GaussianKernel(double x, double bandwidth) 
{
    return std::exp(-0.5 * x * x / (bandwidth * bandwidth)) / (bandwidth * std::sqrt(2.0 * M_PI));
}

double KDE_Estimate(const std::vector<double>& residuals, double t_i, double bandwidth) 
{
    double sum = 0.0;
    for (const double& t_j : residuals) sum += GaussianKernel(t_i - t_j, bandwidth);
    return sum / residuals.size();
}

double ComputeResidualTime(std::vector<double> vertexPos, double originTime, Hit lHit)
{
	int iPMT = lHit.PMT;
	double HitT = (fTimeCorrection + lHit.T) / TimeDelta::ns;
	PMTInfo lPMTInfo = (*fPMTList)[iPMT];
	double distance = Astro_GetDistance(lPMTInfo.Position, vertexPos);

	double tof = distance / fLightSpeed;
	return HitT - tof - originTime;
}

VtxCandidate CreateCandidate(const std::vector<double>& vertex, double lowerLimit, double upperLimit, bool computeSNR)
{
	VtxCandidate out;

	// Copy the vertex info
	out.X = vertex[0];
	out.Y = vertex[1];
	out.Z = vertex[2];
	out.T = vertex[3];
	out.NLL = 0.0; 

	if(computeSNR)
	{
		// We'll do a simple approach: count how many hits are in [lower,upper]
		// for each PMT type, then subtract the dark rate.
	
		int nHitsTotal = fHitCollection->Size();
		int inTimeBnL = 0;
		int inTimemPMT = 0;
	
		for(int iHit = 0; iHit < nHitsTotal; iHit++) {
			const Hit& hit = fHitCollection->At(iHit);
			int iPMT = hit.PMT;
			int pmtType = Astro_GetPMTType(iPMT); // 0=BnL, 1=mPMT
	
			double hitTime = (fTimeCorrection + hit.T) / TimeDelta::ns;
			PMTInfo pInfo = (*fPMTList)[iPMT];
			double distance = Astro_GetDistance(pInfo.Position, vertex);
			double tof = distance / fLightSpeed;
			double residual = hitTime - tof - vertex[3];
			//std::cout << "residual: " << residual << ", " << "lowerlimit and upperlimit " << lowerLimit << ", " << upperLimit << std::endl; 
			if(residual >= lowerLimit && residual <= upperLimit)
			{
				if(pmtType == 0) {
				// BnL
				inTimeBnL++;
				} else {
				// mPMT
				inTimemPMT++;	
				}
			}
		}
	
		double dt = (upperLimit - lowerLimit);
	
		// BnL:
		double darkRateBnL = fTimeWindowSizeFull * fDarkRate_ns[NormalPMT];
		// average dark over dt:
		double darkInWindowBnL = (dt / double(fTimeWindowSizeFull)) * darkRateBnL;
		double signalBnL = std::max(double(inTimeBnL) - darkInWindowBnL, 0.0);
		//std::cout << "inTimeBnL - darkInWindowBnL : " << inTimeBnL << "-" << darkInWindowBnL << std::endl;
		out.SNR = (darkInWindowBnL > 1e-9) ? (signalBnL / darkInWindowBnL) : 999.0;
	}
	return out;
}

void VectorVertexPMT(std::vector<double> vertex, int iPMT, double *dAngles)
{
	// Guillaume 2020/05/20:
	// mPMT referencial is computed once for all PMT in LoadPMTInfo() and MakeMPMTReferencial(iPMT)
	// Guillaume 2020/11/20:
	// Reference is moved to outside LEAF (HKManager or hk-AstroAnalysis)

	PMTInfo lPMTInfo = (*fPMTList)[iPMT];
	int iPMTTop = lPMTInfo.mPMT_RefTube;
	PMTInfo lPMTInfoTop = (*fPMTList)[iPMTTop];

	// 5. Now we have our referential, we should just calculate the angles of the PMT to vertex position vector in this referential.
	// a. calculate the PMT to vertex position vector.

	double dVtx_PMTRef[3];

	dVtx_PMTRef[0] = vertex[0] - lPMTInfoTop.Position[0];
	dVtx_PMTRef[1] = vertex[1] - lPMTInfoTop.Position[1];
	dVtx_PMTRef[2] = vertex[2] - lPMTInfoTop.Position[2];

	double dLengthVtx = Astro_GetLength(dVtx_PMTRef);
	GeoTools::Normalize(dVtx_PMTRef);

	if (VERBOSE >= 3)
	{
		std::cout << "Vertex position = " << vertex[0] << ", " << vertex[1] << ", " << vertex[2] << std::endl;
		std::cout << "Vertex position from PMT = " << dVtx_PMTRef[0] << ", " << dVtx_PMTRef[1] << ", " << dVtx_PMTRef[2] << std::endl;
	}

	// b. Then extract Theta and Phi:
	double dCosTheta = Astro_GetScalarProd(dVtx_PMTRef, lPMTInfo.mPMT_RefZ);
	double dTheta = TMath::ACos(dCosTheta);

	double dPhi = 0.;

	if ((*fPMTList)[iPMT].mPMT_TubeNum == HKAA::kmPMT_TopID)
	{
		// Phi is not defined in that case..
	}
	else
	{
		// We know x=cosPhi x sinTheta and y=sinPhi x sinTheta
		double dX = Astro_GetScalarProd(dVtx_PMTRef, lPMTInfo.mPMT_RefX);
		double dY = Astro_GetScalarProd(dVtx_PMTRef, lPMTInfo.mPMT_RefY);
		double dTanPhi = dY / dX;
		dPhi = TMath::ATan(dTanPhi);

		// tan is symetric from -pi/2 to +pi/2
		if (dX == 0)
		{
			if (dY < 0)
				dPhi = -TMath::Pi() / 2.;
			else
				dPhi = TMath::Pi() / 2.;
		}
		if (dX < 0)
			dPhi += TMath::Pi(); // With this, angle become defines between -pi/2 to -pi/2 +2pi.

		// We wish to bring this from 0 to 2pi:
		if (dPhi < 0)
			dPhi += 2 * TMath::Pi();

		// Actually, we have a symmetry in Phi between [0,pi] and [pi,2pi]. So, we will just define Phi in [0,pi] modulo pi
		if (dPhi > TMath::Pi())
			dPhi = TMath::Pi() - (dPhi - TMath::Pi());
	}

	dAngles[0] = dPhi * 180. / TMath::Pi();
	dAngles[1] = dTheta * 180. / TMath::Pi();
	dAngles[2] = dLengthVtx;
	dAngles[3] = 1; // calculateWeight(dLengthVtx,PMTradius[pmtType],Theta,verbose);

	if (VERBOSE >= 3)
	{
		std::cout << "Angles Phi = " << dAngles[0] << ", Theta = " << dAngles[1] << std::endl;
	}
}

// std::vector<double> ProjectPointToCylinder(std::vector<double> point, double R, double H) 
// {
// 	std::vector<double> dummy;
// 	return dummy;
// }

std::vector<double> ProjectPointToCylinder(std::vector<double> point, double R, double H) 
{
    double x = point[0];
    double y = point[1];
    double z = point[2];
    double halfHeight = H / 2.0;

    // Distance from the z-axis
    double r_xy = std::sqrt(x * x + y * y);

    // --- Project to side surface ---
    double x_side, y_side;
    if (r_xy == 0.0)
	{
        x_side = R;
        y_side = 0.0; // Arbitrary direction
    }
	else 
	{
        x_side = x * R / r_xy;
        y_side = y * R / r_xy;
    }
    double z_side = std::clamp(z, -halfHeight, halfHeight);

    // --- Project to top cap (z = +halfHeight) ---
    double z_top = halfHeight;
    double x_top = x, y_top = y;
    double r_top = std::sqrt(x * x + y * y);
    if (r_top > R && r_top != 0.0) 
	{
        x_top = x * R / r_top;
        y_top = y * R / r_top;
    }

    // --- Project to bottom cap (z = -halfHeight) ---
    double z_bottom = -halfHeight;
    double x_bottom = x, y_bottom = y;
    double r_bottom = std::sqrt(x * x + y * y);
    if (r_bottom > R && r_bottom != 0.0) 
	{
        x_bottom = x * R / r_bottom;
        y_bottom = y * R / r_bottom;
    }

    // --- Compute squared distances to each projection ---
    auto dist2 = [](double dx, double dy, double dz) 
	{
        return dx * dx + dy * dy + dz * dz;
    };

    double d2_side   = dist2(x - x_side,   y - y_side,   z - z_side);
    double d2_top    = dist2(x - x_top,    y - y_top,    z - z_top);
    double d2_bottom = dist2(x - x_bottom, y - y_bottom, z - z_bottom);

    // --- Return closest projection ---
    if (d2_side <= d2_top && d2_side <= d2_bottom) return {x_side, y_side, z_side};
    else if (d2_top <= d2_bottom) return {x_top, y_top, z_top};
    else return {x_bottom, y_bottom, z_bottom};
}