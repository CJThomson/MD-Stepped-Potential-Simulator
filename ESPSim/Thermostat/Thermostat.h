#pragma once

#include <vector>
#include <math.h>

#include "../Random.h"
#include "../Particle.h"
#include "../Sampler.h"
namespace Thermostat
{
  class Thermostat
  {
  public:
    virtual unsigned int runThermostat(std::vector<Particle>& particles,
				       Sampler::Sampler& sampler) = 0;
    virtual double getThermoTime(unsigned long long eventCount) = 0;
    virtual void initialise(const double setT, const Random& rng) = 0;
    virtual bool is_initialised() = 0;
    virtual const char* getType() const = 0; 
  protected:
    Thermostat() {}; //constructor should never be called except by derived classes
  };

  //Andersen Thermostat
  class Andersen: public Thermostat
  {
  public:
  Andersen(double propThermo, bool autoCtrl) :
    proportionEvents(propThermo), autoControl(autoCtrl),
      updateFreq(100), thermofreq(0.0005),
      thermoCount(0), initialised(false) {};
    virtual const char* getType() const { return "Andersen"; } 
     virtual void initialise(const double setT, const Random& rng)
     {
       thermoCount = 0;
       lastUpdate = 0;
       thermofreq = 0.0005;
       setPointTemp = setT;
       RNG = rng;
       initialised = true;
     }
     virtual bool is_initialised () { return initialised; }
     virtual unsigned int runThermostat(std::vector<Particle>& particles,
					Sampler::Sampler& sampler)
     {
       //increase number of thermostat events
       ++thermoCount;
      //generate a particle to collide with
      size_t particleNo = RNG.var_uniformInt(0, particles.size() - 1);
      Particle& p = particles[particleNo];
      Vector3<double> oldv = p.getV();
      double factor = sqrt(setPointTemp);
      //assign values for each component of the particle's velocity from Gaussian
      for(size_t dim (0); dim < 3; ++dim) 
	p.setV()[dim] = RNG.var_normal() * factor;
      sampler.changeKinetic(0.5 * p.getMass() * (p.getV().lengthSqr() - oldv.lengthSqr()));
      return particleNo;
    }
    virtual double getThermoTime(const unsigned long long eventCount)
    {
      if(autoControl)
	{
 	  if(thermoCount > updateFreq)
	    {
	      if(eventCount - lastUpdate != 0) //check divisior is not zero
		{
		  //the proportion of events that are thermostat events divided by target proportion
		  double eventFactor = (double) thermoCount / 
					((eventCount - lastUpdate) * proportionEvents); 
		  thermofreq *= eventFactor;
		  lastUpdate = eventCount;
		  thermoCount = 0;
		}
	    }
	  double t_min = -thermofreq * log(RNG.var_01());
	  return t_min;
	}
      else
	{
	  thermofreq = proportionEvents;
	  double t_min = thermofreq;
	  return t_min;
	}
    }
  private:
    bool initialised;
    double setPointTemp;
    Random RNG;
    bool autoControl;
    unsigned int thermoCount;
    unsigned int updateFreq;
    unsigned long long lastUpdate;
    double thermofreq;
    double proportionEvents;
  };
  //Temperature Rescaling Thermostat
  class Rescale: public Thermostat //needs to be written
  {
    virtual void runThermostat(std::vector<Particle>& particles, double t)
    {
    }
    virtual double getThermoTime(unsigned long long eventCount)
    {
    }
    
  };
}
