#pragma once
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>

#include "Settings.h"
#include "Particle.h"
#include "Lattice.h"
#include "Stepmap.h"
#include "Stepper/include.h"
#include "Stepper/ContPotential.h"

typedef std::vector<Particle>::iterator it_p;

class Simulator
{
 public:
  //Constructor / Deconstructor
  Simulator () {}
  ~Simulator (){ }


  void loadSettings(int, char*[]);
  void initialise();

 private:
  SimSet simSettings;
  SimProp simProperties;
  std::vector<Particle> particles;
  std::vector<std::pair<double, double> > steps;
  Random RNG;
  void zeroMomentum();
};


  
