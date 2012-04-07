#pragma once
#include<math.h> //allows use of mathematical functions
#include<cstdlib> //allows use of random numbers
#include "Vector3.h" //allows use of vector mathematics
class CParticle
{
 public:
  //Constructor
 CParticle(CVector3 location, CVector3 velocity,
	   double radius_in, double mass_in, int ID):
  r(location), v(velocity), radius(radius_in), mass(mass_in), particleNo(ID), r0(location)
  {}

  //vectors on particle
  CVector3 r; //position vector of the particle
  CVector3 r0;
  CVector3 v; //velocity vector of the particle
  CVector3 a; //acceleration vector of the particle

  //particle constants
  double radius; //radius of the particle
  double mass; //mass of the particle
  int particleNo;
  //functions
  double kineticEnergy()//return the specific energy of the particle
  {
    return v.dotProd(v) * mass * 0.5;
  }
};
