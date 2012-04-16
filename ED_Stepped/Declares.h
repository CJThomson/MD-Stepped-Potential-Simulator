#pragma once
//----Structures----
struct Diffusion {
Diffusion(double D, double t):time(t),coDiff(D){ } //Constructor
  double time;
  double coDiff;
};

struct Results
{
  Results() {
    temperature = 0; pressure_d = 0; pressure_c = 0; potential_d =0; potential_c = 0;}
Results(double T, double Pd, double Pc, double Ud, double Uc) :
  temperature(T), pressure_d(Pd), pressure_c(Pc), potential_d(Ud), potential_c(Uc) {}
  const Results& operator += (const Results &r)
  {
    temperature += r.temperature;
    pressure_d += r.pressure_d;
    pressure_c += r.pressure_c;
    potential_d += r.potential_d;
    potential_c += r.potential_c;
    return *this;
  }
  const Results& operator *= (const double &a)
  {
    temperature *= a;
    pressure_d *= a;
    pressure_c *= a;
    potential_d *= a;
    potential_c *= a;
    return *this;
  }
  double temperature;
  double pressure_d;
  double pressure_c;
  double potential_d;
  double potential_c;
};
//----Program Includes----
#include <vector> //allows use of vector (STL structure)
#include <map> //allows use of maps (STL structure)
#include <set> //allows use of sets (STL structure)
#include <utility> //allows use of std::pair
#include <fstream> //allows writing to a file
#include <iostream> //allows writing to console
#include <math.h> //allows use of mathematics
#include <cstdlib> //allows use of random numbers
#include <algorithm> //allows use of sort function
#include <string> //allows use of strings
#include <sstream> //allows use of ostreamstream
#include "Stepper.h" //allows the generation of stepped potentials
#include "Vector3.h" //allows use of vector mathematics
#include "Particle.h" //Event particle Class
#include "Events.h" //Event class
#include "Logger.h" //allows logging
#include "Random.h" //allows use of random numbers
//----TypeDef----
typedef std::vector<CParticle>::iterator it_particle;

//----Function Declarations----
void applyBC(CVector3&);
int calcCell(CVector3);
int calcNewCell(CParticle&);
double calcCellLeave(CParticle&);
double calcCollisionTime(CParticle&, CParticle&, eventTimes::EventType&);
double calcDiff(std::vector<CParticle>&, double);
double calcKinetic(std::vector<CParticle>&);
double calcPressure(double,double);
double calcPotential();
void calcRadDist(std::vector<CParticle> &);
void calcStep(CParticle& ,CParticle&);
double calcSentinalTime(CParticle&);
double calcTemp(std::vector<CParticle>&);
double calcThermoTime(CRandom&, int&);
double calcVelocity(std::vector<CParticle>&, eventTimes&);
void checkCaptureMap(std::vector<CParticle> &);
double continuousP(double);
void continuousRDF(double);
double continuousU();
void correctVelocity(std::vector<CParticle>&);
void initialise(std::vector<CParticle>&, CRandom&);
void initFromFile (std::vector<CParticle>&);
void initSettings(std::vector<std::pair<double, double> >&);
void initSteps();
void generateNeighbourCells(std::vector<std::set<int> >&);
void generateNeighbourList(std::vector<std::set<int> >&,  std::vector<CParticle>&);
void getEvent(CParticle&,
	      std::vector<CParticle>&,
	      std::vector<std::vector<eventTimes> >&,
	      std::vector<eventTimes>&,
	      std::vector<std::set<int> >&,
	      std::vector<std::set<int> >&);
void freeStream(double);
void resetSim();
void runSimulation(std::vector<Results>&, size_t);
void runThermostat(CParticle&, CRandom&, std::vector<eventTimes>&);
void updatePosition(CParticle&);
void zeroMomentum(std::vector<CParticle>&);
