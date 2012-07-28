#pragma once
#include<math.h> //allows use of mathematical functions
#include "Vector3.h" //allows use of vector mathematics
class CParticle
{
 public:
  //Constructors
 CParticle(CVector3<double> location, CVector3<double> velocity,
	   double radius_in, double mass_in, size_t ID):
  r(location), r0(location), v(velocity), radius(radius_in), mass(mass_in), particleNo(ID), collNo(0), updateTime(0)
    {}

  //vectors on particle
  CVector3<double> r; //position vector of the particle
  CVector3<double> r0; //position vector of the particle
  CVector3<double> v; //velocity vector of the particle

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
    return 0.5 * v.lengthSqr() * mass;
  }
  void reset()
  {
    cellNo = 0;
    collNo = 0;
    nextCell = 0;
    updateTime = 0;
  }
};

