#ifndef UNITS_HPP
#define UNITS_HPP

// ENERGY
static constexpr double MeV = 1.;
static constexpr double eV  = MeV * 1E-3;
static constexpr double GeV = MeV * 1E3;
static constexpr double TeV = MeV * 1E6;

// TIME
static constexpr double ns = 1.;
static constexpr double ps = ns * 1E-3;
static constexpr double us = ns * 1E3;
static constexpr double ms = ns * 1E6;
static constexpr double s  = ns * 1E9;

// DISTANCE
static constexpr double mm = 1.;
static constexpr double cm = mm * 1E1;
static constexpr double m  = mm * 1E3;
static constexpr double km = mm * 1E6;

#endif
