
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
  Sampler(unsigned int noParticles, bool countColl, bool gr) :
    numberParticles(noParticles), collCount(countColl), RDF(gr) {};
    void initialise (const std::vector<Particle>&, 
		     const std::vector<std::pair<double, double> >&,
		     const Stepmap&);
    inline void changePotential (const double deltaU) 
    { kineticEnergy -= deltaU; potentialEnergy += deltaU; }
    inline void changeKinetic (const double deltaKE) 
    { kineticEnergy += deltaKE; potentialEnergy -= deltaKE; }
    inline void changeMomentumFlux (double deltaP) { momentumFlux += deltaP;}
    inline double getU() const { return potentialEnergy.current(); }
    inline double getKE() const  { return kineticEnergy; }
    inline double getT() const { return calcTemp(); }
    void eventCount (int, unsigned int,  bool);
    void freeStream(double);

    void sampleRDF(std::vector<Particle>&);
  private:
    //settings
    bool collCount;
    bool RDF;
    
    TAProperty potentialEnergy;
    unsigned int numberParticles;
    double kineticEnergy;
    double momentumFlux;

    TAProperty temperature;
    std::vector<CollisionCount> eventCounts;

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
