#include "Sampler.h"
namespace Sampler
{
  void Sampler::initialise(const std::vector<Particle>& particles, const std::vector<std::pair<double, double> >& steps, const Stepmap& pairStepMap)
  {
    kineticEnergy = 0;
    potentialEnergy = 0;
    for(unsigned int p1 = 0; p1 < particles.size(); ++p1)
      {
	kineticEnergy += particles[p1].kineticEnergy();
	for(unsigned int p2 = p1 + 1; p2 < particles.size(); ++p2)
	  {
	    int stepID(pairStepMap.getStep(p1,p2));
	    potentialEnergy += steps[stepID].second;
	  }
      }

  }
  void Sampler::eventCount ( int type, unsigned int stepNo, bool inwards )
  {
    if(collCount)
      {
	if(!inwards) { ++stepNo; } //if an outward event increase stepNo to get correct step
	eventCounts[stepNo].incrCount(type, inwards); //increment step Count
      }
  }

  void Sampler::freeStream(double deltaTime)
  {
    temperature = calcTemp();
    temperature.stream(deltaTime);
    potentialEnergy.stream(deltaTime);
  }

  void Sampler::sampleRDF(std::vector<Particle>& particles )
  {
    
    /*for(it_particle p1 = particles.begin(); p1 != particles.end(); ++p1) // loop through every particle
      {
	updatePosition(*p1); //update the position of the particle
	for(it_particle p2 = p1 + 1; p2 != particles.end(); ++p2) //loop through every other particle
	  {
	    updatePosition(*p2); //update the postion of that particle too
	    PBCVector<double> distance(length, true, p1->r - p2->r); //calculate the distance between the particles
	    //applyBC(distance); //apply PBC to the distance
	    if(distance.length() < maxR) //if the distance is less than the max radius considered
	      {
		int index = floor(distance.length() * noBins/ maxR); //calculate which bin the particle is in
		if(index < 0 || index >= noBins) //if an invalid bin then send error
		  {
		    cerr << "ERROR: Invalid index is calcRadDist: " << index << endl;
		    exit(1);
		  }
		++rdf_d[index]; //add 1 to the bin counter
	      }
	  }
	  }*/
  }
}
