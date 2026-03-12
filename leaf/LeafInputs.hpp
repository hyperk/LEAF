#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <mutex>

//* Root Headers
#include "TRandom3.h"

//* WCSim Headers
#include "WCSimRootGeom.hh"
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafConfig.hpp"

extern const Geometry* fGeometry;
extern const HitCollection<Hit>* fHitCollection;
extern TimeDelta fTriggerTime;
extern TimeDelta fTimeCorrection;

extern const std::vector<PMTInfo> *fPMTList;
extern double fDarkRate_ns[NPMT_CONFIGURATION];

extern std::vector< std::vector<double> > fPositionList;

extern int fThread;

extern TRandom3 * fRand;

extern std::vector<double> fTrueVtxPos;
extern std::vector<double> fTrueDir;

void InitInputs(const Geometry *lGeometry);

void MakePositionList();

void SetNThread(int iThread=N_THREAD);

//* True Vertex isn't used for the fit, just to get feedback on how well the fit is at different steps
void SetTrueVertexInfo(std::vector<double> vtx, double time);
void SetTrueDirInfo(std::vector<double> trueDir);


