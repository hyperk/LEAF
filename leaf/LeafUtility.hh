#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

//* DataModel informations
#include "Geometry.h"
#include "HitCollection.h"
#include "LeafConfig.hh"

inline bool SortOutputVector ( const std::vector<double>& v1, const std::vector<double>& v2 ) 
{ 
	return v1[4] < v2[4]; 
}

double Distance3D(std::vector<double> point1, std::vector<double> point2);
		
bool CorrectCandidate(std::vector<std::vector<double>>* candidates, std::vector<double> point, std::vector<double> fTrueVtxPos);

bool ContainsTrueVtx(std::vector<std::vector<double>>* tRecoVtxPos, std::vector<double> fTrueVtxPos);