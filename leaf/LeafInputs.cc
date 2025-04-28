#include "LeafInputs.hh"

const Geometry* fGeometry;

const HitCollection<Hit>* fHitCollection;
TimeDelta fTriggerTime;
TimeDelta fTimeCorrection;

const std::vector<PMTInfo> *fPMTList;
double fDarkRate_ns[NPMT_CONFIGURATION];

void InitInputs(const Geometry *lGeometry)
{
    fGeometry = lGeometry;
    fPMTList = fGeometry->GetPMTList();

    fDarkRate_ns[NormalPMT] = fGeometry->pmt_dark_rate[HKAA::kIDPMT_BnL];
	fDarkRate_ns[MiniPMT] = fGeometry->pmt_dark_rate[HKAA::kIDPMT_3inch];
}