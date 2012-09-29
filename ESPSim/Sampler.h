
#pragma once

#include <vector>
#include <utility>
#include <boost/foreach.hpp>

#include "Particle.h"
#include "TAProperty.h"
#include "Stepmap.h"
namespace Sampler
{
  class CollisionCount;
  
  class Sampler
  {
  public:
  Sampler(unsigned int noParticles, double m, double rho, bool countColl, bool gr) :
    numberParticles(noParticles), mass(m), density(rho), 
      collCount(countColl), RDF(gr) {};
    void initialise (const std::vector<Particle>&, 
		     const std::vector<std::pair<double, double> >&,
		     const Stepmap&);
    inline void changePotential (const double deltaU) 
    { kineticEnergy -= deltaU; potentialEnergy += deltaU; }
    inline void changeKinetic (const double deltaKE) 
    { kineticEnergy += deltaKE; potentialEnergy -= deltaKE; }
    inline void changeMomentumFlux (double deltaP) { momentumFlux += deltaP;}
    inline double getU() const { return potentialEnergy.current() / numberParticles; }
    inline double getKE() const  { return kineticEnergy; }
    inline double getT() const { return calcTemp(); }
    inline double getMeanT() const { return temperature.mean(); }
    inline double getMeanU() const { return potentialEnergy.mean() / numberParticles; }
    inline double getMomFlux() const { return momentumFlux.current(); }
    inline double getFluxTime() const { return momentumFlux.time(); }
    inline double getP() const { return calcPressure(); }
    void eventCount (int, unsigned int,  bool);
    void freeStream(double);

    void sampleRDF(std::vector<Particle>&);
  private:
    //settings
    bool collCount;
    bool RDF;

    unsigned int numberParticles;    
    double mass;
    double density;

    TAProperty potentialEnergy;
    double kineticEnergy;
    TAProperty momentumFlux;
    TAProperty temperature;
    std::vector<CollisionCount> eventCounts;
    
    inline double calcPressure() const
    { 
      double pressure = density * temperature.mean() + mass * density * momentumFlux.current() 
	/ (numberParticles * 3.0 * momentumFlux.time()); 

      if(abs(pressure) > 1E10)
	{
	  std::cerr << density << " - " << temperature.mean() << " - " << mass
		    << " - " << density << " - " << momentumFlux.mean() 
		    << " - " << momentumFlux.current()
		    << " - " << momentumFlux.time();
	  exit(0);
	}
      return pressure;
    }


    inline double calcTemp() const
    { return kineticEnergy / (1.5 * numberParticles); }
  };

class CollisionCount
  {
  public:
    CollisionCount() { reset(); }
    void reset()
    {
      _Core = 0;
      _Capture = 0;
      _Release = 0;
      _Disassociation = 0;
      _Association = 0;
      _BounceIn = 0;
      _BounceOut = 0;
    }
    void incrCount(int type, bool inwards)
    {
      switch(type)
	{
	case 1: //Core event
	  ++_Core;
	  break;
	case 2: 
	  if(inwards) //Capture event
	    ++_Capture;
	  else //Release event
	    ++_Release;
	    break;
	case 3:
	  if(inwards) //Association
	    ++_Association;
	  else //Dissociation
	    ++_Disassociation;
	    break;
	case 4: //bounce (in/out)
	  if(inwards)
	    ++_BounceIn;
	  else
	    ++_BounceOut;
	  break;
	default: //error
	  break;
	}
    }

  private:
    unsigned int _Core;
    unsigned int _Capture;
    unsigned int _Release;
    unsigned int _Disassociation;
    unsigned int _Association;
    unsigned int _BounceIn;
    unsigned int _BounceOut;
  };
}
