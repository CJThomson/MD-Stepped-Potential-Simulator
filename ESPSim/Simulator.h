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

  inline const std::vector<Particle>& getParticles() const { return particles; }
  inline const Stepmap& getStepMap() const { return stepmap; }
  inline const std::vector<std::pair<double, double> >& getSteps() const { return steps; }
  inline double getSysLength() { return simProperties.getLength(); }
  bool isRunning(double, unsigned long long , bool);
  double progress(double, unsigned long long , bool);

  inline std::vector<Particle>& setParticles() { return particles; }
  inline Stepmap& setStepMap() { return stepmap; }
  inline SimSet& setSettings() { return simSettings; }
  inline SimProp& setProperties() { return simProperties; }
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


  
