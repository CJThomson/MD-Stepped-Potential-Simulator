#pragma once

#include <math.h>

class ContPotential
{
 public:
  virtual double energy(double r) = 0;
  virtual double force(double r) = 0;
 protected:
  ContPotential() {};
};

class LennardJones: public ContPotential
{
 public:
  LennardJones(double s, double e) :
   lj_sig(s), lj_eps(e) {};
  virtual inline double energy(double r)
  {
    return 4.0 * lj_eps * (pow(lj_sig / r, 12) - pow(lj_sig / r, 6));
  }
  virtual inline double force(double r)
  {
    return 24.0 * lj_eps / lj_sig * (2.0 * pow(lj_sig / r, 13) - pow(lj_sig / r, 7));
  }
 private: 
  double lj_eps;
  double lj_sig;
};
