#pragma once

#include <vector>
#include <utility>
#include <math.h>
#include <boost/shared_ptr.hpp>

#include "Stepper.h"
#include "ContPotential.h"
#include "../Maths/NumTech.h"
#include "../Maths/Functions.h"
#include "../Maths/polyRoot.h"

typedef std::vector<std::pair<double, double> >::iterator it_step;

namespace Stepper
{
  class Enr_Virial: public Stepper
  {
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
  public:
  Enr_Virial(double temp, boost::shared_ptr<ContPotential> pot) :
    Stepper(0, 0, temp, 0, 0, pot) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      Functions* func = new Partition(lambda, potential);
      {
	Maths::NumTech integrator;
	for(it_step i = steps.begin(); i != steps.end(); ++i)
	  {	  
	    double energy;
	    if(i != steps.end() - 1)
	      {
		it_step j = i + 1;
		energy = - lambda 
		  * log(3.0 / ( 4.0 * M_PI * (pow(i->first, 3) - pow(j->first, 3)))
			* integrator.integrator(func, j->first, i->first, 100));
	      }
	    else
	      {
		energy = - lambda 
		  * log(3.0 / ( 4.0 * M_PI * (pow(i->first, 3)))
			* integrator.integrator(func, 1e-16, i->first, 100));
	      }
	    i->second = energy;
	  }
	delete func;
      }
    }
  };
  class Enr_Mid: public Stepper
  {
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
  public:
  Enr_Mid(boost::shared_ptr<ContPotential> pot) :
    Stepper(0, 0, 0, 0, 0, pot) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    {
	      it_step j = i + 1;
	      energy = potential->energy((i->first + j->first) * 0.5);
	    }
	  else
	    {
	      energy = potential->energy(i->first);
	    }
	  i->second = energy;
	}
    }
  };
  class Enr_Average: public Stepper
  {
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
  public:
  Enr_Average(boost::shared_ptr<ContPotential> pot) :
    Stepper(0, 0, 0, 0, 0, pot) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      Functions* func = new Partition(lambda, potential);
      {
	Maths::NumTech integrator;
	for(it_step i = steps.begin(); i != steps.end(); ++i)
	  {	  
	    double energy;
	    if(i != steps.end() - 1)
	      {
		it_step j = i + 1;
		double v1 = pow(i->first,3);
		double v2 = pow(j->first,3);
		double vtotal = 1 / (v1 + v2);
		energy = v1 * vtotal * potential->energy(i->first)
		  + v2 * vtotal * potential->energy(j->first);
	      }
	    else
	      {
		energy = potential->energy(i->first);
	      }
	    i->second = energy;
	  }
	delete func;
      }
    }
  };
}
