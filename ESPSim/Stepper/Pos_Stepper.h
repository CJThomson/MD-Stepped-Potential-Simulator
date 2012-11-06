#pragma once

#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>

#include "Stepper.h"
#include "ContPotential.h"
#include "../Maths/NumTech.h"
#include "../Maths/Functor.h"
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
    double rCut;
    unsigned int noSteps;
  public:
  Pos_Even(double _rCut, unsigned int _noSteps) : rCut(_rCut), noSteps(_noSteps){};
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
    double T;
    double rCut;
    double noSteps;
    boost::shared_ptr<ContPotential> potential;
  public:
  Pos_ExpForce(double _T, double _rCut, unsigned int _noSteps, 
	       boost::shared_ptr<ContPotential> _potential) :
    T(_T), rCut(_rCut), noSteps(_noSteps), potential(_potential) {};

    virtual void genPositions(std::vector<std::pair<double, double> >& steps)
    {
      boost::shared_ptr<Functor> f(new ExpForce(T, potential) );
      Maths::NumTech integrator;
      steps.clear();
      steps.resize(noSteps);
      double totalEF = integrator.integrator(f, ZERO, rCut, 100);
      double r_lower = ZERO;
      for(size_t i(1); i <= noSteps; ++i)
	{
	  boost::shared_ptr<Functor> f2( new IntegralEq(f, 100, r_lower) );
	  double step = integrator.rootFinder(f2, r_lower, rCut, 
					      1e3, 1e-5, totalEF / noSteps);
	  if(step != 0)
	    steps[noSteps - i].first = r_lower = step;
	}
    }
 
  };
  class Pos_EvenEnergy: public Stepper
  {
  public:
  Pos_EvenEnergy(double _energyInt, double _rCut, 
		 boost::shared_ptr<ContPotential> _potential) :
    energyInt(_energyInt), rCut(_rCut), potential(_potential) {}

    virtual void genPositions(std::vector<std::pair<double, double> >& steps)
    {
      Polynomial rootFinder;
      double maxEnergy = 0;
      size_t i(0);
      if(energyInt == 0)
	{
	  std::cerr << "Error: Energy interval of 0 passed to EvenEnergy Stepper" << std::endl;
	  exit(3);
	}
      //define polynomial of Lennard Jones
      std::vector<double> polynomial(13);
      polynomial[6] = -4;
      polynomial[12] = 4;

      std::vector<double> roots;
      steps.clear();
      steps.push_back(std::make_pair(rCut,0)); //add the cutoff radius
      while(maxEnergy < 50)
	{
	  polynomial[0] = - (-1 + 0.5 * energyInt  + energyInt * i) - potential->shift();
	  maxEnergy = -polynomial[0];
	  rootFinder.rootFind(polynomial, roots, 1e-5, 1e2);
	  for(std::vector<double>::iterator j = roots.begin(); j != roots.end(); ++j)
	    if(*j >= 0 && *j < rCut)
	      if(find(steps.begin(), steps.end(), 
		      std::pair<double, double>(*j, 0))==steps.end())
		steps.push_back(std::make_pair(*j, 0));
	  if(maxEnergy <= 0)
	    sort(steps.begin(), steps.end(), std::greater<std::pair<double, double> >());
	  ++i;
	}
    }
  private:
    bool sortPair(std::pair<double, double> a, std::pair<double, double> b)
    { return a.first > b.first; }
    double energyInt;
    double rCut;
    boost::shared_ptr<ContPotential> potential;
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps) {}
  };


}
