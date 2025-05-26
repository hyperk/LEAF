#include "LeafInputs.hh"

const Geometry* fGeometry;

const HitCollection<Hit>* fHitCollection;
TimeDelta fTriggerTime;
TimeDelta fTimeCorrection;

const std::vector<PMTInfo> *fPMTList;
double fDarkRate_ns[NPMT_CONFIGURATION];

int fThread = N_THREAD;

TRandom3 * fRand;

std::vector<double> fTrueVtxPos;
std::vector<double> fTrueDir;

std::vector< std::vector<double> > fPositionList;

void InitInputs(const Geometry *lGeometry)
{
    fGeometry = lGeometry;
    fPMTList = fGeometry->GetPMTList();
    fThread = N_THREAD;
    fRand = new TRandom3();

    fDarkRate_ns[NormalPMT] = fGeometry->pmt_dark_rate[HKAA::kIDPMT_BnL];
	fDarkRate_ns[MiniPMT] = fGeometry->pmt_dark_rate[HKAA::kIDPMT_3inch];
    MakePositionList();
}

void MakePositionList()
{
	// Make position list for a given step size
	double dStep = fSearchVtxStep;
    fPositionList = std::vector< std::vector<double> >();
	fPositionList.clear();

	for (double time = -50; time < 50; time += (dStep / fLightSpeed))
	{
		for (double radius = 0; radius <= fTankRadius; radius += dStep)
		{
			double perimeter = 2 * TMath::Pi() * radius;
			int numberOfStepsOnCircle = floor(perimeter / dStep);
			if (numberOfStepsOnCircle == 0) numberOfStepsOnCircle = 1;
			double angleStepOnCircle = 2 * TMath::Pi() / numberOfStepsOnCircle;

			for (double angle = 0; angle <= 2 * TMath::Pi(); angle += angleStepOnCircle)
			{
				for (double height = -fTankHalfHeight; height <= fTankHalfHeight; height += dStep)
				{
					std::vector<double> tVtxPosition(4, 0.); // In centimeters
					tVtxPosition[0] = radius * TMath::Cos(angle);
					tVtxPosition[1] = radius * TMath::Sin(angle);
					tVtxPosition[2] = height;
					tVtxPosition[3] = time;
					fPositionList.push_back(tVtxPosition);
				}
			}
		}
	}
}

void SetTrueVertexInfo(std::vector<double> vtx, double time) 
{
	fTrueVtxPos = std::vector<double>(5, 0.);
	fTrueVtxPos[0] = vtx[0];
	fTrueVtxPos[1] = vtx[1];
	fTrueVtxPos[2] = vtx[2];
	fTrueVtxPos[3] = time;
	fTrueVtxPos[4] = 0.;
}

void SetTrueDirInfo(std::vector<double> trueDir)
{
	fTrueDir = std::vector<double>(3, 0.);
	fTrueDir[0] = trueDir[0];
	fTrueDir[1] = trueDir[1];
	fTrueDir[2] = trueDir[2];
}

void SetNThread(int iThread){ fThread=iThread; }