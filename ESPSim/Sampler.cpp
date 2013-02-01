#include "Sampler.h"
namespace Sampler
{
  void Sampler::initialise(const std::vector<Particle>& particles, const std::vector<std::pair<double, double> >& steps, const Stepmap& pairStepMap)
  {
    kineticEnergy = 0;
    potentialEnergy = 0;
    temperature = 0;
    momentumFlux = 0;
    for(unsigned int p1 = 0; p1 < particles.size(); ++p1)
      {
	kineticEnergy += particles[p1].kineticEnergy();
	for(unsigned int p2 = p1 + 1; p2 < particles.size(); ++p2)
	  {
	    int stepID(pairStepMap.getStep(p1,p2)); 
	    if(stepID != -1)
	      potentialEnergy += steps[stepID].second;
	  }
      }

    if(collCount)
      for(size_t i(0); i < steps.size(); ++i )
	eventCounts.push_back(CollisionCount(i, steps[i].first, steps[i].second));
    time(&startTime);
  }
  void Sampler::terminate(const unsigned long long eventCount, const double t)
  {
    time(&endTime);

    double timeDiff = difftime(endTime, startTime);
    eventPS = eventCount / timeDiff;
    timePS = t / timeDiff;
    eventCnt = eventCount;
    simTime = t;
  }

  void Sampler::eventCount (const int type, int stepNo, const bool inwards )
  {
    if(collCount)
      {
	if(inwards) { ++stepNo; } //if an inward event increase stepNo to get correct step
	eventCounts[stepNo].incrCount(type, inwards); //increment step Count
      }
  }

  void Sampler::freeStream(double deltaTime)
  {
    temperature = calcTemp();
    temperature.stream(deltaTime);
    potentialEnergy.stream(deltaTime);
    momentumFlux.stream(deltaTime);
  }
  void Sampler::initialiseRDF(const double noBins, const double maxR, const double timeInt)
  {
    RDF_timeInt = timeInt;
    RDF_maxR = maxR;
    RDF_bins = noBins;
    RDF_noReadings = 0;
    RDF_data.clear();
    RDF_data.resize(noBins);
  }
  
  void Sampler::sampleRDF(const std::vector<Particle>& particles, const double sysLength )
  {
    ++RDF_noReadings;
    double bin_per_width = RDF_bins / RDF_maxR;
    for(unsigned int p1 = 0; p1 < particles.size(); ++p1)
      for(unsigned int p2 = p1 + 1; p2 < particles.size(); ++p2)
	{
	  //Calculate the distance between two particles
	  PBCVector<double> distance(sysLength, true, particles[p1].getR() - particles[p2].getR());
	  if(distance.length() < RDF_maxR)
	    {
	      //Calculate which bin the particle pair is in
	      int index = floor(distance.length() * bin_per_width);
	      //and increment value
	      ++RDF_data[index];
	    }
	}
  }
  
  std::vector<double> Sampler::calcRDF(unsigned int numberParticles, double density) const
  {
    double deltaR = RDF_maxR / RDF_bins;
    std::vector<double> RDF_return;
    RDF_return.resize(RDF_data.size());
    for(size_t i = 0; i < RDF_bins; ++ i)
      {
	double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
	RDF_return[i] = RDF_data[i] / (0.5 * numberParticles * RDF_noReadings * volShell * density);  
      }
    return RDF_return;
  }
}
