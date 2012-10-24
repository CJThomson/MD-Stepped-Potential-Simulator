#pragma once

#include <math.h>

class ContPotential
{
 public:
  virtual double energy(double r) = 0;
  virtual double force(double r) = 0;
  virtual double shift() = 0;
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
  virtual inline double shift()
  {
    return 0;
  }
 private: 
  double lj_eps;
  double lj_sig;
};

class LennardJonesShifted: public ContPotential
{
 public:
 LennardJonesShifted(double s, double e, double rCut) :
  lj_sig(s), lj_eps(e), r_cut(rCut)
  { lj_shift = 4.0 * lj_eps * (pow(lj_sig / r_cut, 12) - pow(lj_sig / r_cut, 6)); }
  virtual inline double energy(double r)
  {
    return 4.0 * lj_eps * (pow(lj_sig / r, 12) - pow(lj_sig / r, 6))
      - lj_shift;
  }
  virtual inline double force(double r)
  {
    return 24.0 * lj_eps / lj_sig * (2.0 * pow(lj_sig / r, 13) - pow(lj_sig / r, 7));
  }
  virtual inline double shift()
  {
    return lj_shift;
  }

 private:
  double lj_eps;
  double lj_sig;
  double r_cut;
  double lj_shift;
};
