
#pragma once

#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include <ctime>

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
    void terminate (const unsigned long long, const double);
    inline void changePotential (const double deltaU) 
    { potentialEnergy += deltaU; }
    inline void changeKinetic (const double deltaKE) 
    { kineticEnergy += deltaKE; }
    inline void changeMomentumFlux (double deltaP) { momentumFlux += deltaP;}
    inline double getU() const { return potentialEnergy.current() / numberParticles; }
    inline double getKE() const  { return kineticEnergy; }
    inline double getT() const { return calcTemp(); }
    inline double getMeanT() const { return temperature.mean(); }
    inline double getMeanSqrT() const {return temperature.meanSqr(); }
    inline double getMeanU() const { return potentialEnergy.mean() / numberParticles; }
    inline double getMeanSqrU() const 
    {return potentialEnergy.meanSqr() / (numberParticles * numberParticles); }
    inline double getMomFlux() const { return momentumFlux.current(); }
    inline double getFluxTime() const { return momentumFlux.time(); }
    inline double getP() const { return calcPressure(); }
    inline double getEventPS() const { return eventPS; }
    inline double getTimePS() const { return timePS; }
    inline char* getStartTime() const 
    { char* retVal = asctime(localtime(&startTime));
      retVal[strlen(retVal)-1] = '\0';
      return retVal; }
    inline char* getEndTime() const 
    { char* retVal = asctime(localtime(&endTime));
      retVal[strlen(retVal)-1] = '\0';
      return retVal; }
    inline double getTime() const { return simTime; }
    inline unsigned int getEventCount() const { return eventCnt; }
    inline double getMeanFreeTime() const { return simTime / eventCnt; }
    inline const std::vector<CollisionCount>& getCollCount() const { return eventCounts; }
    inline const bool getRDF() const { return RDF; }
    void eventCount (const int, int,  const bool);
    void freeStream(double);
    void initialiseRDF(const double noBins, const double maxR, const double timeInt);
    void sampleRDF(const std::vector<Particle>&, const double);
    std::vector<double> calcRDF(const unsigned int, const double) const;
    inline double getBinWidth() const { return RDF_maxR / RDF_bins; }
    inline double getRDFTime(const double time) const { return time + RDF_timeInt; }
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
    double eventPS;
    double timePS;
    time_t startTime;
    time_t endTime;
    unsigned int eventCnt;
    double simTime;
    
    double RDF_timeInt;
    double RDF_maxR;
    unsigned int RDF_bins;
    unsigned int RDF_noReadings;
    
    std::vector<unsigned int> RDF_data;
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
  CollisionCount(unsigned int id, double pos, double enr) :
    _ID(id), _r(pos), _U(enr)
      {
	reset();
      }
    void reset()
    {
      _Capture = 0;
      _Release = 0;
      _BounceIn = 0;
      _BounceOut = 0;
    }
    unsigned int getID() const { return _ID; }
    double getR() const { return _r; }
    double getU() const { return _U; }
    unsigned int getCount(int type, bool inwards) const
    {
      switch(type)
	{
	case 2: 
	  if(inwards) //Capture event
	    return _Capture;
	  else //Release event
	    return _Release;
	    break;
	case 4: //bounce (in/out)
	  if(inwards)
	    return _BounceIn;
	  else
	    return _BounceOut;
	  break;
	default: //error
	  break;
	}
    }
    void incrCount(int type, bool inwards)
    {
      switch(type)
	{
	case 2: 
	  if(inwards) //Capture event
	    ++_Capture;
	  else //Release event
	    ++_Release;
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
    unsigned int _ID;
    double _r;
    double _U;
    unsigned int _Capture;
    unsigned int _Release;
    unsigned int _BounceIn;
    unsigned int _BounceOut;
  };
}
