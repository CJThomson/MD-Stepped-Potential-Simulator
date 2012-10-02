#pragma once

#include <math.h>
#include <boost/shared_ptr.hpp>
#include "../Stepper/ContPotential.h" 

class Functions
{
 public:
  virtual double f(double r) = 0;
 protected:
 Functions(double T, boost::shared_ptr<ContPotential> pot) :
  beta(1/T), potential(pot) {};
  double beta;
  boost::shared_ptr<ContPotential> potential;

};

class Partition: public Functions
{
 public:
 Partition(double T, boost::shared_ptr<ContPotential> pot) :
  Functions(T, pot) {};

  virtual inline double f(double r)
  {
    return 4 * M_PI * exp(-beta * potential->energy(r)) * r * r;
  }
};

class ExpForce: public Functions
{
 public:
 ExpForce(double T, boost::shared_ptr<ContPotential> pot) :
  Functions(T, pot) {};

  virtual inline double f(double r)
  {
    return fabs(potential->force(r) * 4 * M_PI * exp(-beta * potential->energy(r)) * r * r);
  }
};

