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
  bool sampleRDF;
  unsigned int RDF_bins;
  double RDF_maxR;
  double RDF_timeInt;
  unsigned long long simEvents;
  unsigned long long eqEvents;
  unsigned int noRuns;
  double simTime;
  double eqTime;
  bool thermoOn;
  bool thermoControl;
  double thermoFreq;
  bool reduceOut;
  boost::shared_ptr<Thermostat::Thermostat> thermostat;
  std::string nlType;
 public:
  //constructor
 SimSet(): reduceOut(0){}

  //get access 
  inline int getOutLog() const { return writeOutLog; }
  inline unsigned long long getRunEvent() const { return simEvents; }
  inline unsigned long long getEQEvent() const { return eqEvents; }
  inline double getRunTime() const { return simTime; }
  inline double getEQTime() const { return eqTime; }
  inline boost::shared_ptr<Thermostat::Thermostat> getThermostat()
  { return thermostat; }
  inline bool activeThermo() const { return thermoOn; }
  inline bool getThermoControl() const {return thermoControl; }
  inline double getThermoFreq() const { return thermoFreq; }
  inline bool getSampleColl() const { return sampleColl; }
  inline bool getSampleRDF() const { return sampleRDF; }
  inline unsigned int getRDF_bins() const { return RDF_bins; }
  inline double getRDF_maxR() const { return RDF_maxR; }
  inline double getRDF_timeInt() const { return RDF_timeInt; }
  inline unsigned int getRuns() const { return noRuns; }
  inline bool getReducedOut() const { return reduceOut; }
  inline const char* getNLType() const { return nlType.c_str(); }
  bool isTime(bool eq) { return (eq) ? eqTime != 0 : simTime != 0;  }
  //set access
  void setOutLog(int value) { writeOutLog = value; }
  void setSampleColl(bool value) { sampleColl = value; }
  void setSampleRDF(unsigned int bins, double maxR, double timeInt) 
  { 
    sampleRDF = true; RDF_bins = bins; RDF_maxR = maxR; RDF_timeInt = timeInt;
  }

  void setRunEvent(unsigned long long value) { simEvents = value; }
  void setEQEvent(unsigned long long value) { eqEvents = value; }
  void setRunTime(double value) { simTime = value; }
  void setEQTime(double value) { eqTime = value; }
  void setRuns(unsigned int value) { noRuns = value; }
  void setThermoControl(bool value) { thermoControl = value; }
  void setThermoFreq(double value) { thermoFreq = value;  }
  void setNLType(const char* value) { nlType = value; }
  void setReducedOut(bool value) { reduceOut = value; }
  void setThermoType(const char* type) 
  { 
    if(strcmp(type, "None") == 0)
      {
	thermoOn = false;
    	boost::shared_ptr<Thermostat::Thermostat> 
	  tempThermo(new Thermostat::None()); 
	thermostat = tempThermo;
      }
    else if(strcmp(type,"Andersen") == 0 )
      { 
	thermoOn = true;
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
      ("noRuns,R", boost::program_options::value<unsigned int>(),
       "Number of Production Runs to Simulate")
      ("NL", boost::program_options::value<int>(),
       "Type of neighbour list: \n\t0: None\n\t1: Simple")
      ("ReducedOut", boost::program_options::value<bool>(),
       "Suppress the output freq for cluster scripting")
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
    if(vm.count("noRuns"))
      noRuns = vm["noRuns"].as<unsigned int>();
    if(vm.count("ReducedOut"))
      reduceOut = vm["ReducedOut"].as<bool>();
    if(vm.count("NL"))
      { switch (vm["NL"].as<int>()) { 
	case 0: nlType = "None"; break;
	case 1: nlType = "Simple"; break; }
      }
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
       "Stepped potential positions: \n\t0: Even              1: Even Energy\n\t2: Exp Mean Force    3: Chapela")
      ("stepenr", boost::program_options::value<int>(),
       "Stepped potential energies: \n\t0: Mid     1: Left     2: Right\n\t3: Average Volume      4: Volume Averaged Energy\n\t5: Expected Energy     6: Virial     7: Chapela")
      ("stepcore", boost::program_options::value<int>(),
       "Stepped potential core: \n\t0: None\n\t1: Manual")
      ("corepos", boost::program_options::value<double>(),
       "Position of core in stepped potential")
      ("nosteps,n", boost::program_options::value<unsigned int>(),
       "Number of discontinuities in stepped potential")
      ("energyInt,u", boost::program_options::value<double>(),
       "Energy interval between steps (EvenEnergy only)")
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
	case 1: stepPositions = "EvenEnergy"; break; 
	case 3: stepPositions = "Chapela"; break; }
      }
    if(vm.count("stepenr"))
      { switch (vm["stepenr"].as<int>()) { 
	case 0: stepEnergies = "Mid"; break;
	case 1: stepEnergies = "Left"; break;
	case 2: stepEnergies = "Right"; break;
	case 3: stepEnergies = "AvgVol"; break; 
	case 4: stepEnergies = "AvgEnr"; break; 
	case 5: stepEnergies = "ExpEnr"; break; 
	case 6: stepEnergies = "Virial"; break; 
	case 7: stepEnergies = "Chapela"; break; }
      }
    if(vm.count("nosteps"))
      noSteps = vm["nosteps"].as<unsigned int>();

    if(vm.count("energyInt"))
      energyInterval = vm["energyInt"].as<double>();
    if(vm.count("energyInt"))
      energyInterval = vm["energyInt"].as<double>();
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
