#pragma once

#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>

#include "Stepper.h"
#include "ContPotential.h"
#include "../Maths/NumTech.h"
#include "../Maths/Functions.h"
#include "../Maths/polyRoot.h"

#define ZERO 1e-16
namespace Stepper
{
  class Pos_Even: public Stepper
  {
  private:
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps) {}
  public:
  Pos_Even(double rCut, unsigned int n) :
    Stepper(0,0, 0, rCut, n) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      for(size_t i = noSteps; i > 0; --i)
	steps.push_back(std::make_pair(i * rCut / noSteps, 0));
    } 

  };


  class Pos_ExpForce: public Stepper
  {
  private:
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps) {}
  public:
  Pos_ExpForce(double temp, double rCut, unsigned int n, boost::shared_ptr<ContPotential> pot) :
    Stepper(0, 0, temp, rCut, n, pot) {};

    virtual void genPositions(std::vector<std::pair<double, double> >& steps)
    {
      Functions* func = new ExpForce(lambda, potential);
      {
	Maths::NumTech integrator;
	steps.clear();
	steps.resize(noSteps);
	double totalEF = integrator.integrator(func, ZERO, rCut, 100);
	double r_lower = ZERO;
	for(size_t i(1); i <= noSteps; ++i)
	  {
	    double step = integrator.limitSolver(func, totalEF / noSteps,
						 r_lower, rCut, 100, 1e3, 1e-5);
	    if(step != 0)
	      steps[noSteps - i].first = r_lower = step;
	  }
	delete func;
      }
    }
 
  };
  class Pos_EvenEnergy: public Stepper
  {
  public:
  Pos_EvenEnergy(double energyStep, double rCut, boost::shared_ptr<ContPotential> pot) :
    Stepper(0,0, energyStep, rCut, 0, pot) {};

    virtual void genPositions(std::vector<std::pair<double, double> >& steps)
    {
      Polynomial rootFinder;
      double maxEnergy;
      size_t i(2);
      std::vector<double> polynomial(13);
      polynomial[6] = -4;
      polynomial[12] = 4;
      std::vector<double> roots;

      steps.clear();
      steps.resize(noSteps);  
      steps.begin()->first = rCut; //add the cut off radius
      while(maxEnergy < 200)
	{
	  polynomial[0] = - (-1 + 0.5 * lambda  + lambda * i);
	  maxEnergy = polynomial[0];
	  rootFinder.rootFind(polynomial, roots, 1e-5, 1e2);
	  for(std::vector<double>::iterator j = roots.begin(); j != roots.end(); ++j)
	    if(*j >= 0 && *j < rCut)
	      if(find(steps.begin(), steps.end(), 
		      std::pair<double, double>(*j, 0))==steps.end())
		steps[noSteps - i].first = *j;
	  if(maxEnergy <= 0)
	    sort(steps.begin(), steps.end(), std::greater<std::pair<double, double> >());
	  ++i;
	}
    }
  private:
    bool sortPair(std::pair<double, double> a, std::pair<double, double> b)
    { return a.first > b.first; }
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps) {}
  };


}
