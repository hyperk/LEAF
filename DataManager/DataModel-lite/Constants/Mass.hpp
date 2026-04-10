#ifndef MASS_HPP
#define MASS_HPP

#include <map>

#include "Constants/Units.hpp"

static const std::map<int, float> PDGParticleMasses {
	{  11, 0.510998910 * MeV}, //   electron
	{2212,  938.272013 * MeV}, //   proton
	{2112,   939.56536 * MeV}  //   neutron
};
#endif
