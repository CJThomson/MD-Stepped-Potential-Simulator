#pragma once
#include <vector>
#include "Vector3.h"
#include "Particle.h"

class Lattice
{
 public:
  virtual void placeParticles(std::vector<Particle>& particles, double systemLen) = 0;
 protected:
  Lattice(){} // no copies of this class should be created.

};


class FCC:public Lattice
{
 public:
  virtual void placeParticles(std::vector<Particle>& particles, double length)
  {
    int n = ceil(pow(particles.size()/ 4.0, 1.0 / 3.0)); //find cubic root of number of particles
    double a = length / n; //calculate the scale factor
    int j = 0, x = 0, y = 0, z = 0;
    for(size_t i (0); i < particles.size(); ++i) //for every particle
      {
	Vector3<double> location;
	switch (j) //sort the particles between the different particles in the FCC root vector
	  {
	  case 0:
	    location = Vector3<double>(x * a - length * 0.5, y * a - length * 0.5, z * a - length * 0.5);
	    break;
	  case 1:
	    location = Vector3<double>(x* a - length * 0.5, (y + 0.5) * a - length * 0.5, (z + 0.5) * a - length * 0.5);
	    break;
	  case 2:
	    location = Vector3<double>((x + 0.5) * a - length * 0.5, y * a - length * 0.5, (z + 0.5) * a - length * 0.5);
	    break;
	  case 3:
	    location = Vector3<double>((x + 0.5) * a - length * 0.5, (y + 0.5) * a - length * 0.5, z * a - length * 0.5);
	    j = -1;
	    break;
	  }
	++j; 
	if(j == 0) //once a particle has been placed in every spot in the unit vector
	  ++x; //more to the right
	if(x >= n)
	  {
	    x=0; //then move up
	    ++y;
	  }
	if(y >= n) //then move forward
	  {
	    y=0;
	    ++z;
	  }

	particles[i].setR() = location;
	particles[i].setID(i);
      }

  }
};
