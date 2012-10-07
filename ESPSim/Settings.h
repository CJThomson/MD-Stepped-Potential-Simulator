#pragma once
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <math.h>
#include <string>

#include "Thermostat/Thermostat.h"
class Settings
{
 public:
  virtual void setOptions(boost::program_options::options_description& simOpts) = 0;
  virtual void loadCLSettings(boost::program_options::variables_map& vm) = 0;

 protected:
  void optionWarning (const char* message)
  { std::cerr << "Warning: " << message << std::endl; }
  void optionError (const char* message)
  { std::cerr << "Error (2): " << message << " Sim will now close." << std::endl; exit(2);}
};

class SimSet: public Settings
{
 private:
  int writeOutLog;
  bool runForTime;
  bool sampleColl;
  unsigned long long simEvents;
  unsigned long long eqEvents;
  double simTime;
  double eqTime;
  bool thermoControl;
  double thermoFreq;
  boost::shared_ptr<Thermostat::Thermostat> thermostat;
 public:
  //constructor
 SimSet(): writeOutLog(-1){}

  //get access 
  inline int getOutLog() const { return writeOutLog; }
  inline unsigned long long getRunEvent() const { return simEvents; }
  inline unsigned long long getEQEvent() const { return eqEvents; }
  inline double getRunTime() const { return simTime; }
  inline double getEQTime() const { return eqTime; }
  inline boost::shared_ptr<Thermostat::Thermostat> getThermostat()
  { return thermostat; }
  inline bool getThermoControl() const {return thermoControl; }
  inline double getThermoFreq() const { return thermoFreq; }
  inline bool getSampleColl() const { return sampleColl; }
  bool isTime(bool eq) { return (eq) ? eqTime != 0 : simTime != 0;  }
  //set access
  void setOutLog(int value) { writeOutLog = value; }
  void setSampleColl(bool value) { sampleColl = value; }
  void setRunEvent(unsigned long long value) { simEvents = value; }
  void setEQEvent(unsigned long long value) { eqEvents = value; }
  void setRunTime(double value) { simTime = value; }
  void setEQTime(double value) { eqTime = value; }
  void setThermoControl(bool value) { thermoControl = value; }
  void setThermoFreq(double value) { thermoFreq = value;  }
  void setThermoType(const char* type) 
  { 
    if(strcmp(type,"Andersen") == 0 )
      { 
	boost::shared_ptr<Thermostat::Thermostat> 
	  tempThermo(new Thermostat::Andersen(thermoFreq, thermoControl)); 
	thermostat = tempThermo;
      }
    else
      { std::cerr << type<< "- Error: Invalid thermostat type"; exit(4); }
 }
  virtual void setOptions(boost::program_options::options_description& simOpts)
  {
    simOpts.add_options()
      ("eqevents", boost::program_options::value<unsigned int>(),
       "Equilibration Length (events)")
      ("runevents", boost::program_options::value<unsigned int>(),
       "Production Run Length (events)")
      ("eqtime", boost::program_options::value<double>(),
       "Equilibration Length (time)")
      ("runtime", boost::program_options::value<double>(),
       "Production Run Length (time)")
      /*("thermotype", boost::program_options::value<int>(),
	"Thermostat Type: \n\t1: Andersen" << - Add at some point*/
      ;
  }

  virtual void loadCLSettings(boost::program_options::variables_map& vm)
  {
    if(vm.count("eqevents"))
      eqEvents = vm["eqevents"].as<unsigned int>();
    if(vm.count("eqtime"))
      eqTime = vm["eqtime"].as<double>();
    if(vm.count("runevents"))
      simEvents = vm["runevents"].as<unsigned int>();
    if(vm.count("runtime"))
      simTime = vm["runtime"].as<double>();


  }

};
class SimPotential: public Settings
{
 public:
  SimPotential() {};
  virtual void setOptions(boost::program_options::options_description& simOpts)
  {
    simOpts.add_options()
      ("potential", boost::program_options::value<int>(),
       "Continuous potential: \n\t0: Lennard Jones\n\t1: Shifted LJ")
      ("epsilon,e", boost::program_options::value<double>(),
       "Potential minimum")
      ("sigma,s", boost::program_options::value<double>(),
       "Potential root")
      ("rcut,r", boost::program_options::value<double>(),
       "Cut-off radius")
      ("steppos", boost::program_options::value<int>(),
       "Stepped potential positions: \n\t0: Even\n\t1: Even Energy\n\t2: Exp Mean Force")
      ("stepenr", boost::program_options::value<int>(),
       "Stepped potential energies: \n\t0: Mid Values\n\t1: Virial")
      ("stepcore", boost::program_options::value<int>(),
       "Stepped potential core: \n\t0: None\n\t1: Manual")
      ("corepos", boost::program_options::value<double>(),
       "Position of core in stepped potential")
      ("nosteps,n", boost::program_options::value<unsigned int>(),
       "Number of discontinuities in stepped potential")
      ("energyinterval,u", boost::program_options::value<double>(),
       "Energy interval between steps (EvenEnergy only")
      ;
  }

  virtual void loadCLSettings(boost::program_options::variables_map& vm)
  {
    if(vm.count("potential"))
      { switch (vm["potential"].as<int>()) {
	case 0: contPotential = "LennardJones"; break;
	case 1: contPotential = "LennardJones_shifted"; break; }
      }
    if(vm.count("epsilon"))
      epsilon = vm["epsilon"].as<double>();
    if(vm.count ("sigma"))
      epsilon = vm["sigma"].as<double>();
    if(vm.count("rcut"))
      rCutOff = vm["rcut"].as<double>();
    if(vm.count("steppos"))
      { switch (vm["steppos"].as<int>()) { 
	case 0: stepPositions = "Even"; break;
	case 1: stepPositions = "EvenEnergy"; break; }
      }
    if(vm.count("stepenr"))
      { switch (vm["stepenr"].as<int>()) { 
	case 0: stepEnergies = "Mid"; break;
	case 1: stepEnergies = "Virial"; break;
	case 2: stepEnergies = "Average"; break; }
      }
    if(vm.count("nosteps"))
      noSteps = vm["nosteps"].as<unsigned int>();
    if(vm.count("energyinterval"))
      noSteps = vm["energyInterval"].as<double>();
    
  }
  //get
  inline const char* getContPotential() const { return contPotential.c_str(); }
  inline const char* getStepPositions() const { return stepPositions.c_str(); }
  inline const char* getStepEnergies() const { return stepEnergies.c_str(); }
  inline const char* getStepCore() const { return stepCore.c_str(); }
  inline const char* getStepPotential() const { return stepPotential.c_str(); }
  inline const unsigned int getNoStep() const { return noSteps; }
  inline const double getEpsilon() const { return epsilon; }
  inline const double getSigma() const { return sigma; }
  inline const double getRCutOff() const { return rCutOff; }
  inline const double getEnergyInterval() const { return energyInterval; }
  inline const double getCorePos() const { return corePos; }
  inline const int getCoreNoSigm() const { return coreNoSigma; }

  //set
  void setContPotential(const char* type, const double cut, const double e, const double r)
  {
    contPotential = type;
    rCutOff = cut;
    epsilon = e;
    sigma = r;
  }
  void setStepPositions(const char* type, const unsigned int n) 
  { stepPositions = type; noSteps = n; }
  void setStepPositions(const char* type, const double de) 
  { stepPositions = type; energyInterval = de; }
  void setStepEnergies(const char* type) 
  { stepEnergies = type; }
  void setStepCore(const char* type, double pos = 0) 
  { stepCore = type; corePos = pos; }
  void setStepCore(const char* type, unsigned int sigma) 
  { stepCore = type; coreNoSigma = sigma; }
  void setStepPotential(const char* type)
  { stepPotential = type; }

 private:
  std::string contPotential;
  std::string stepPositions;
  std::string stepEnergies;
  std::string stepCore;
  std::string stepPotential;
  unsigned int noSteps;
  double epsilon;
  double sigma;
  double rCutOff;
  double energyInterval;
  double corePos;
  double coreNoSigma;
  
};
class SimProp:public Settings
{
 private:
  double temperature;
  double density;
  unsigned int numberOfParticles;
  double systemLength;

  
 public:
  //constructor - just use default
 SimProp() : temperature (0), density(0), numberOfParticles(0), systemLength(0)
    {}

  //get access
  double getT() const { return temperature; }
  double getDensity() const { return density; }
  unsigned int getN() const { return numberOfParticles; }
  double getLength() 
  { 
    //check if systemLength has been calculated before
    if(systemLength == 0) { systemLength = pow(numberOfParticles / density, 1.0 / 3.0);}
    return systemLength;
  }
  void setT(double value) { temperature = value; }
  void setDensity(double value) { density = value; }
  void setN(unsigned int value) { numberOfParticles = value; }
  virtual void setOptions(boost::program_options::options_description& simOpts)
  {
    namespace po = boost::program_options;
    simOpts.add_options()
      ("temperature,T", po::value<double>(), "System temperature")
      ("density,d", po::value<double>(), "System Density")
      ("particles,N", po::value<unsigned int>(), "Number of particles")
      ;
  }

  virtual void loadCLSettings(boost::program_options::variables_map& vm)
  {
    if(vm.count("temperature"))
      temperature = vm["temperature"].as<double>();

    if(temperature == 0)
      optionError("No temperature specified");

    if(vm.count("density"))
      density = vm["density"].as<double>();
    
    if(density == 0)
      optionError("No density specified");

    if(vm.count("particles"))
      numberOfParticles = vm["particles"].as<unsigned int>();
    if(numberOfParticles == 0)
      optionError("Number of particles not specified");
  }

};
