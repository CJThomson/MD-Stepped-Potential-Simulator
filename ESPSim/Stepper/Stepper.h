//This is the new and improved stepper...It will be awesome and OO
#pragma once
#include <vector>
#include <utility>
#include <limits>
#include <boost/shared_ptr.hpp>

#include "ContPotential.h"
#include "../Maths/Functor.h"
namespace Stepper
{
  class Stepper
  {
  public:
    virtual void addCore(std::vector<std::pair<double, double> >& steps) = 0;
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) = 0;
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) = 0;
    virtual void genPotential(std::vector<std::pair<double, double> >& steps) = 0;
  protected:
    Stepper() {};
  Stepper(double s, double e, double l, double rCutoff, unsigned int n, 
	  boost::shared_ptr<ContPotential> pot = boost::shared_ptr<ContPotential> ()) : 
    sigma(s), epsilon(e), lambda(l),rCut(rCutoff), noSteps(n), potential(pot) 
    {};

    double sigma; //holds particle diameter
    double epsilon; //holds well energy
    double lambda; //holds misc other values, length of well (square well) energy interval (energy stepping), temperature
    double rCut; //holds the cut off radius
    boost::shared_ptr<ContPotential> potential; //the continuous potential to be used
    unsigned int noSteps;
  };

}

