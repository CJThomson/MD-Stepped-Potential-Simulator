#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include "Settings.h"
#include "Particle.h"
#include "Lattice.h"
#include "Logger.h"
#include "Stepmap.h"
#include "ParseXML.h"
#include "Stepper/include.h"
#include "Stepper/ContPotential.h"
#include "Thermostat/Thermostat.h"

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
  void productionRun();

  inline const std::vector<Particle>& getParticles() const { return particles; }
  inline const Stepmap& getStepMap() const { return stepmap; }
  inline const std::vector<std::pair<double, double> >& getSteps() const { return steps; }
  inline double getSysLength() { return simProperties.getLength(); }
  inline double getDensity() const { return simProperties.getDensity(); } 
  inline double getTemperature() const { return simProperties.getT(); }
  inline const SimSet& getSettings() { return simSettings; }
  inline const Random& getRNG() const { return RNG; }
  bool isRunning(double, unsigned long long , bool);
  double progress(double, unsigned long long , bool);

  inline std::vector<Particle>& setParticles() { return particles; }
  inline boost::shared_ptr<Thermostat::Thermostat> setThermostat()
  { return simSettings.getThermostat(); }  
  inline Stepmap& setStepMap() { return stepmap; }
  inline SimSet& setSettings() { return simSettings; }
  inline SimProp& setProperties() { return simProperties; }
 private:
  SimSet simSettings;
  SimProp simProperties;
  SimPotential simPot;
  std::vector<Particle> particles;
  std::vector<std::pair<double, double> > steps;
  Random RNG;
  Stepmap stepmap;

  void zeroMomentum();
  bool runSim();
  void initSteps();

};


  
