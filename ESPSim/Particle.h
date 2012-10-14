#pragma once
#include "Vector3.h" //allows use of vector mathematics
#include "Random.h"
class Particle
{

 private:
  Vector3<double> r; //position vector of the particle
  Vector3<double> r0; //position vector of the particle
  Vector3<double> v; //velocity vector of the particle
  //particle constants
  double radius; //radius of the particle
  double mass; //mass of the particle
  unsigned int id;
  unsigned int cellNo;
  unsigned long long noColl;
  unsigned int nextCell;
  double lastUpdate;
 public:
  //Constructors
 Particle(Vector3<double> location, Vector3<double> velocity,
	  double radius_in, double mass_in, unsigned int particleID):
  r(location), r0(location), v(velocity), radius(radius_in), mass(mass_in), id(particleID), noColl(0), lastUpdate(0), cellNo(0)
    {}
 Particle() : r(), r0(), v(), radius(1), mass(1), id(-1), noColl(0), lastUpdate(0){ }
  //GET access functions
  inline const Vector3<double>& getR() const { return r; }
  inline const Vector3<double>& getR0() const { return r0; }
  inline const Vector3<double>& getV() const { return v; }

  inline double getRadius() const { return radius; }
  inline double getMass() const { return mass; }
  inline unsigned int getID() const { return id; }
  unsigned int getCell() const { return cellNo; }
  unsigned long long getNoColl() const { return noColl; }
  unsigned int getNextCell() const { return nextCell; }
  double getLastUpdate() const { return lastUpdate; }

  //SET access functions
  void setR(Vector3<double> &newr) { r = newr;}
  void setV(Vector3<double> &newv) { v = newv;}
  Vector3<double>& setR() { return r; }
  Vector3<double>& setV() { return v; }
  void resetR0() {r0 = r;}
  void setID(unsigned int value) { id = value; }
  void setNextCell(int newCell) {nextCell = newCell;}
  void setLastUpdate(double currentTime) {lastUpdate = currentTime;}

  void incrNoColl() { ++noColl; }

  //functions
  double kineticEnergy() const //return the specific energy of the particle
  {
    return 0.5 * v.lengthSqr() * mass;
  }
  void reset()
  {
    cellNo = 0;
    noColl = 0;
    nextCell = 0;
    lastUpdate = 0;
  }

  void move(const double sysTime)
  {
    std::cerr << "still being called" << std::endl;
    if(sysTime != lastUpdate)
      {
	r += v * (sysTime - lastUpdate);
	lastUpdate = sysTime;
      }
  }
  void resetV(Random& rng)
  {
    for(size_t dim(0); dim < 3; ++dim)
      v[dim] = rng.var_normal(); //give each particle a normally distributed velocity
  }
};

