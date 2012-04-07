#pragma once
//header with all includes, function and type declarations for the ForceSim Program
//----Structures----
struct Diffusion {
  Diffusion(double D, double t):time(t),coDiff(D){ } //Constructor
  double time;
  double coDiff;
};

//----Includes----
#include <vector> //allows use of vector (storage structure)
#include <fstream> //allows writing to file
#include <iostream> //allows writing to screen
#include <math.h> //allows use of mathematics
#include <cstdlib> //allows use of 'random' numbers
#include <iomanip> //allows use of setprecision
//include boost libraries for normal distributions
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>

#include "Vector3.h" //allows use of vector mathematics
#include "Particle.h" //allows use of the particle class
#include "Logger.h" //allows logging of files

//----Function Declarations----
void applyBC(CVector3&);
void calcAllForces(std::vector<CParticle> &,std::vector<int> &, int []);
double calcDiff(std::vector<CParticle>&, double);
CVector3 calcForce(CVector3);
double calcKinetic(std::vector<CParticle> &);
double calcPotential(std::vector<CParticle> &, std::vector<int>&, int[]);
double calcVirial(std::vector<CParticle> &, std::vector<int>&, int[]);
void calcRadDist(std::vector<CParticle> &);
double calcTemp(std::vector<CParticle>&);
void correctVelocity(std::vector<CParticle>&);
void initialise(std::vector<CParticle>&);
void initFromFile(std::vector<CParticle>&);
double calcBoltzmannH(std::vector<CParticle> &);
void zeroMomentum(std::vector<CParticle> &);
void calcNeighbourList(std::vector<CParticle> &, std::vector<int> &, int []);
//Define pi
const double PI = 3.1415926536;
