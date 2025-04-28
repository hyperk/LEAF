#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>

//WCSim Headers
#include "WCSimRootGeom.hh"
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafConfig.hh"

extern const Geometry* fGeometry;
extern const HitCollection<Hit>* fHitCollection;
extern TimeDelta fTriggerTime;
extern TimeDelta fTimeCorrection;

extern const std::vector<PMTInfo> *fPMTList;
extern double fDarkRate_ns[NPMT_CONFIGURATION];

    // 1.385;//1.373;//refraction index of water

void InitInputs(const Geometry *lGeometry);


