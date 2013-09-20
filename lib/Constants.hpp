//Sean M. Law

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <map>
#include <string>

const double PI = 4.0*atan(1.0);
const double DEG2RAD = PI/180.0;
const double RAD2DEG = 180.0/PI;

//Reference: CRC Handbook for Chemistry and Physics, 1983/1984
//CHARMM Source: source/ltm/consta_ltm.src
//SI Units

const double ONEKCAL = 4184.0; // 1Kcal = 4184 J
const double NAVO = 6.022045E23; //Avogadro constant, 1/mol
const double KBOLTZ = (1.380662E-23); //Boltzmann constant, J/K
const double AMU = 1.6605655E-27; //Atomic mass unit, kg
const double RGAS = 8.314472; //Molar gas constant, J/K
const double PLANCK = 6.6260693E-34; //Planck constant, J*s

//AKMA Units
const double AKMATPS = 4.88882129E-02; //One AKMA time unit as picoseconds
const double AKMATS = 4.88882129E-14; //One AKMA time unit as seconds
const double kB = KBOLTZ*NAVO/ONEKCAL; //Boltzmann constant, kcal/mol/K
const double Rgas = 8.314472*NAVO/ONEKCAL; //Molar gas constant kcal/mol/K
const double planck = PLANCK*NAVO/ONEKCAL; //Planck constant, (kcal/mol)*s


#endif
