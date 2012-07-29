#include <iostream>
#include <boost/program_options.hpp>
#include <math.h>

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

 public:
  //constructor
 SimSet(): writeOutLog(-1){}

  //get access
  int getOutLOg() { return writeOutLog; }

  //set access
  void setOutLog(int value) { writeOutLog = value; }

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
 SimProp() : temperature (0), density(0), numberOfParticles(0), systemLength(0){}

  //get access
  double getT() { return temperature; }
  double getDensity() { return density; }
  unsigned int getN() { return numberOfParticles; }
  double getLength() 
  { 
    //check if systemLength has been calculated before
    if(systemLength == 0) { pow(numberOfParticles / density, 1.0 / 3.0);}
    return systemLength;
  }

  virtual void setOptions(boost::program_options::options_description& simOpts)
  {
    namespace po = boost::program_options;
    simOpts.add_options()
      ("temperature,T", po::value<double>(), "System temperature")
      ("density, rho", po::value<double>(), "System Density")
      ("particles, N", po::value<unsigned int>(), "Number of particles")
      ;
  }

  virtual void loadCLSettings(boost::program_options::variables_map& vm)
  {
    if(vm.count("temperature") && temperature == 0)
      temperature = vm["temperature"].as<double>();
    else
      optionError("No temperature specified");

    if(vm.count("density") && density == 0)
      density = vm["density"].as<double>();
    else
      optionError("No density specified");

    if(vm.count("particles") && numberOfParticles == 0)
      numberOfParticles = vm["particles"].as<unsigned int>();
    else
      optionError("Number of particles not specified");
  }

};
