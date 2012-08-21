#pragma once
#include <map>
#include <utility>
#include <algorithm>
#include <math.h>

#include "Vector3.h"
#include "Particle.h"

typedef std::vector<Particle>::iterator it_p;

class Stepmap
{
 public:
  void populateMap(std::vector<Particle>& particles, 
		   std::vector<std::pair<double, double> >& steps, double sysLen)
  {
    for(it_p p1 = particles.begin(); p1 != particles.end(); ++p1)
      for(it_p p2 = p1 + 1; p2 != particles.end(); ++p2)
	calcStep(*p1, *p2, steps, sysLen);
  }  

  bool checkMap(std::vector<Particle>& particles, 
		std::vector<std::pair<double, double> >& steps, double sysLen)
  {
    bool check = checkCaptureMap(particles, steps, sysLen);
    if(!check)
      populateMap(particles, steps, sysLen);
    return check;
  }
  void deleteFromMap()
  {
     pairStepMap.erase(it_last); 
  }
  void moveInwards(unsigned int p1, unsigned int p2)
  {
    unsigned int i = std::min(p1, p2);
    unsigned int j = std::max(p1, p2);
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>::iterator it_map 
      = pairStepMap.find(std::make_pair(i, j)); //find collision state of particles
    if(it_map == pairStepMap.end()) 
      pairStepMap.insert(std::make_pair(std::make_pair(p1, p2), 0)); 
    else
      ++(it_map->second); //move particles in one step
  }
  void moveOutwards(unsigned int p1, unsigned int p2)
  {
    unsigned int i = std::min(p1, p2);
    unsigned int j = std::max(p1, p2);
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>::iterator it_map 
      = pairStepMap.find(std::make_pair(i, j)); //find collision state of particles
    if(it_map->second == 0)
      pairStepMap.erase(it_map);
    else
      --(it_map->second);
  }
  int getStep(unsigned int p1, unsigned int p2)
  {
    unsigned int i = std::min(p1, p2);
    unsigned int j = std::max(p1, p2);
    std::map<std::pair<unsigned int, unsigned int>, unsigned int>::iterator it_map 
      = pairStepMap.find(std::make_pair(i, j)); //find collision state of particles
    it_last = it_map;
    if(it_map == pairStepMap.end())
      return -1;
    return it_map->second;
  }
 private:
  std::map<std::pair<unsigned int, unsigned int>, unsigned int> pairStepMap;
  std::map<std::pair<unsigned int, unsigned int>, unsigned int>::iterator it_last;
  void calcStep(Particle& particle1, Particle& particle2, 
		std::vector<std::pair<double, double> >& steps, double length)
  {
    using namespace std;
    unsigned int p1 = particle1.getID();
    unsigned int p2 = particle2.getID();
    PBCVector<double> r12(length, true, particle1.getR() - particle2.getR()); 
    double distance = r12.length(); //calculate the length of the separation vector
    unsigned int size = steps.size();
    //if distance is outside the map leave function
    if(distance > steps[0].first) 
      return;
    for(size_t i(1); i < size; ++i) //loop through all the possible steps
      {
	if(distance > steps[i].first) //if the distance is in the step
	  {
	    //insert into the step map
	    pairStepMap.insert(make_pair(make_pair(p1, p2), i - 1)); 
	    break; //break the loop
	  }
	else if(distance < steps[size - 1].first)
	  {
	    //insert into the step map
	    pairStepMap.insert(make_pair(make_pair(p1, p2), size - 1)); 
	    break; //break the loop
	  }	
      }
  }

  bool checkCaptureMap(std::vector<Particle> &particles,
		       std::vector<std::pair<double, double> >& steps, double length)
  {
    //CAPTURE TEST

    for(int i = 0; i < particles.size();++i)
      for(int j = i + 1; j < particles.size(); ++j)
	{
	  PBCVector<double> r12(length, true, particles[i].getR()- particles[j].getR());
	  std::map<std::pair<unsigned int, unsigned int>, unsigned int>::const_iterator it_map 
	    = pairStepMap.find(std::pair<int,int>(i, j)); //find collision state of particles
	  double distance = r12.length();
	  int size = steps.size();
	  bool invalid = false;

	  //if particles aren't pairstepmap but should be
	  if(it_map == pairStepMap.end() && distance <= steps[0].first)
	    invalid = true;

	  //if particles are in pairstepmap but shouldn't be
	  if(it_map != pairStepMap.end() && distance > steps[0].first)
	    invalid = true;
	  if(it_map != pairStepMap.end()) //make sure step is in capture map

	    {
	      //if particles are in a step further out than they should be
	      if (steps[it_map->second].first < distance)
		invalid = true;
	      
	      //if paritcles are in a step further in than they should be
	      if(it_map->second != size - 1)
		if(steps[it_map->second + 1].first > distance)
		  invalid = true;
	    }
	  if(invalid)
	    {
	      std::cerr << "\rWarning: Capture map invalid. Regenerating...\n" ;
	      return false;
	    }
	}
    return true;
  }
};
  
