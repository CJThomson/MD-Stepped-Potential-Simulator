#pragma once
#include<math.h> //allows use of mathematical functions
#include "Vector3.h" //allows use of vector mathematics
class CParticle
{
 public:
  //Constructors
 CParticle(CVector3 location, CVector3 velocity,
	   double radius_in, double mass_in, int ID):
  r(location), r0(location), v(velocity), radius(radius_in), mass(mass_in), particleNo(ID), collNo(0), updateTime(0)
    {}

  //vectors on particle
  CVector3 r; //position vector of the particle
  CVector3 r0; //position vector of the particle
  CVector3 v; //velocity vector of the particle

  //particle constants
  double radius; //radius of the particle
  double mass; //mass of the particle
  int particleNo;
  int cellNo;
  int collNo;
  int nextCell;
  double updateTime;
  //functions
  double kineticEnergy()//return the specific energy of the particle
  {
    return 0.5*v.dotProd(v)*mass;
  }
};

