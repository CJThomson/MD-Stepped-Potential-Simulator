#pragma once
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>

#include "Settings.h"
#include "Particle.h"
#include "Lattice.h"
#include "Stepmap.h"
#include "ParseXML.h"
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
  void equilibrate();

  std::vector<Particle> getParticles() const { return particles; }
  Stepmap getStepMap() const { return stepmap; }
  std::vector<std::pair<double, double> > getSteps() const { return steps; }
  double getSysLength() { return simProperties.getLength(); }
  bool isRunning(double, unsigned long long , bool);
  double progress(double, unsigned long long , bool);

  std::vector<Particle>& setParticles() { return particles; }
  Stepmap& setStepMap() { return stepmap; }
  SimSet setSettings() { return simSettings; }
  SimProp setProperties() { return simProperties; }
 private:
  SimSet simSettings;
  SimProp simProperties;
  std::vector<Particle> particles;
  std::vector<std::pair<double, double> > steps;
  Random RNG;
  Stepmap stepmap;

  void zeroMomentum();
  bool runSim();

};


  
