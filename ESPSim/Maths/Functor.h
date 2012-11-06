#pragma once

#include <math.h>
#include <boost/shared_ptr.hpp>
#include "../Stepper/ContPotential.h" 


class Functor
{
 public:
  virtual double operator()(double r) = 0;
 protected:
  Functor() {}
};

class Partition: public Functor
{
 public:
  Partition(double _T, boost::shared_ptr<ContPotential> _potential) :
  T(_T), potential(_potential) {};

  virtual inline double operator()(double r)
  {
    double beta = 1 / T;
    return 4 * M_PI * exp(-beta * potential->energy(r)) * r * r;
  }
 private:
  double T;
  boost::shared_ptr<ContPotential> potential;
};

class ExpForce: public Functor
{
 public:
 ExpForce(double _T, boost::shared_ptr<ContPotential> _potential) :
  T(_T), potential(_potential) {};

  virtual inline double operator()(double r)
  {
    double beta = 1 /T;
    return fabs(potential->force(r) * 4 * M_PI * exp(-beta * potential->energy(r)) * r * r);
  }
 private:
  double T;
  boost::shared_ptr<ContPotential> potential;
};
#include "NumTech.h"
class IntegralEq: public Functor
{
 public:
 IntegralEq(boost::shared_ptr<Functor> _f, unsigned int _intervals,
	    double _lower) :
  f(_f), intervals(_intervals), lower(_lower) {}
  virtual inline double operator()(double r)
  {
    Maths::NumTech integrator;
    return integrator.integrator(f, lower, r, intervals);
  }
 private:
  boost::shared_ptr<Functor> f;
  unsigned int intervals;
  double lower;
};

