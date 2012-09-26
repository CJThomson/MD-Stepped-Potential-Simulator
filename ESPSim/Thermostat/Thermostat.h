#pragma once

#include <vector>
#include <math.h>

#include "Random.h"
#include "Particle.h"
namespace Thermostat
{
  class Thermostat
  {
  public:
    virtual void runThermostat(std::vector<Particle>& particles, double t, unsigned int p1=-1) = 0;
    virtual double getThermoTime(bool autoControl) = 0;
    virtual int getParticle() = 0;
  protected:
    Thermostat() {}; //constructor should never be called except by derived classes
  };

  //Andersen Thermostat
  class Andersen: public Thermostat
  {
  public:
  Andersen(double setT, Random& rng) :
    setPointTemp(setT), RNG(rng) {};
    virtual void runThermostat(std::vector<Particle>& particles, double t, unsigned int p1=-1)
    {
      ++thermoCount;
      CVector3<double> oldv = p.getV();
      p.move(t);

      double factor = sqrt(setPointTemp);
      //assign values for each component of the particle's velocity from Gaussian
      for(size_t dim (0); dim < 3; ++dim) 
	particle.v[dim] = RNG.var_normal() * factor;
    }
    virtual double getThermoTime(bool autoControl)
    {

      if(thermoCount > thermoUpdate)
	{
	  if((eventCount - thermoLastUpdate != 0)) //check divisior is not zero
	    {
	      thermoMeanFreeTime *= (double) thermoCount 
		/ ((eventCount - thermoLastUpdate) * thermoSetting);
	      thermoLastUpdate = eventCount;
	      thermoCount = 0;
	    }
	}

      double t_min = -thermoMeanFreeTime * log(RNG.var_01());
      return t_min;
    }
  private:
    double setPointTemp;
    Random RNG;
    unsigned int thermoCount;
  };
  //Temperature Rescaling Thermostat
  class Rescale: public Thermostat
  {
    
  }
}
