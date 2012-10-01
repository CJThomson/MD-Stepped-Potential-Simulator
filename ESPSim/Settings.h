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
  unsigned long long simEvents;
  unsigned long long eqEvents;
  double simTime;
  double eqTime;
  bool thermoControl;
  bool thermoFreq;
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
  bool isTime(bool eq) { return (eq) ? eqTime != 0 : simTime != 0;  }
  //set access
  void setOutLog(int value) { writeOutLog = value; }
  void setRunEvent(unsigned long long value) { simEvents = value; }
  void setEQEvent(unsigned long long value) { eqEvents = value; }
  void setRunTime(double value) { simTime = value; }
  void setEQTime(double value) { eqTime = value; }
  void setThermoControl(bool value) { thermoControl = value; }
  void setThermoFreq(double value) { thermoFreq = value; }
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
      ("outlog", boost::program_options::value<int>()->default_value(0), 
       "level of output logging: \n\t0:no logging\n\t1:event descriptions\n\t2:full descriptions")
      ;
  }

  virtual void loadCLSettings(boost::program_options::variables_map& vm)
  {
    if(vm.count("outlog") && writeOutLog == -1)
      writeOutLog = vm["outlog"].as<int>();
    else
      optionWarning("no output logging level specified, default value of 0 used.");
  }

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
