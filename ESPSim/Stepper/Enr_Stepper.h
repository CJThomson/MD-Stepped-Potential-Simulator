#pragma once

#include <vector>
#include <utility>
#include <math.h>
#include <boost/shared_ptr.hpp>

#include "Stepper.h"
#include "ContPotential.h"
#include "../Maths/NumTech.h"
#include "../Maths/Functor.h"
#include "../Maths/polyRoot.h"

typedef std::vector<std::pair<double, double> >::iterator it_step;

namespace Stepper
{
  class Enr_Mid: public Stepper
  {
  public:
  Enr_Mid(boost::shared_ptr<ContPotential> _potential) :
    potential(_potential) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    energy = potential->energy((i->first + (i + 1)->first) * 0.5);
	  else
	    energy = potential->energy(i->first);
	  i->second = energy;
	}
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_Left: public Stepper
  {
  public:
  Enr_Left(boost::shared_ptr<ContPotential> _potential) :
    potential(_potential) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    energy = potential->energy((i + 1)->first);
	  else
	    energy = potential->energy(i->first);
	  i->second = energy;
	}
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_Right: public Stepper
  {
  public:
  Enr_Right(boost::shared_ptr<ContPotential> _potential) :
    potential(_potential) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	i->second = potential->energy(i->first);
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_Virial: public Stepper
  {
  public:
  Enr_Virial(double _T, boost::shared_ptr<ContPotential> _potential) :
    T(_T), potential(_potential) {};
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      boost::shared_ptr<Functor> func (new Partition(T, potential) );
      Maths::NumTech integrator;
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    {
	      it_step j = i + 1;
	      energy = -T 
		* log(3.0 / ( 4.0 * M_PI * (pow(i->first, 3) - pow(j->first, 3)))
		      * integrator.integrator(func, j->first, i->first, 100));
	    }
	  else
	    {
	      energy = - T 
		* log(3.0 / ( 4.0 * M_PI * (pow(i->first, 3)))
		      * integrator.integrator(func, 1e-16, i->first, 100));
	    }
	  i->second = energy;
	}
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    double T;
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_AvgVol: public Stepper
  {

  public:
  Enr_AvgVol(boost::shared_ptr<ContPotential> _potential) :
    potential(_potential) {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    {
	      double ravg = pow(0.5 * ( pow(i->first, 3) + pow((i + 1)->first, 3) ), 
				1.0 / 3.0 );
	      energy = energy = potential->energy(ravg);
	    }
	  else
	    energy = potential->energy(i->first);
	  i->second = energy;
	}
      
    }
  private:
    boost::shared_ptr<ContPotential> potential;
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
  };
  class Enr_AvgEnr: public Stepper
  {
  public:
  Enr_AvgEnr(boost::shared_ptr<ContPotential> _potential) :
    potential(_potential) {};
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      boost::shared_ptr<Functor> func (new PotentialVolInt(potential) );
      Maths::NumTech integrator;
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    {
	      double volume = 4.0 / 3.0 * M_PI *
		(pow(i->first, 3) - pow((i + 1)->first, 3));
	      energy = integrator.integrator(func, 
					     (i + 1)->first, i->first,
					     100);
	      energy /= volume;
	    }
	  else
	    {
	      double volume = 4.0 / 3.0 * M_PI * pow(i->first, 3);
	      energy = integrator.integrator(func, 
					     0, i->first,
					     100);
	      energy /= volume;
	    }
	  i->second = energy;
	}
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_ExpEnr: public Stepper
  {
  public:
  Enr_ExpEnr(double _T, boost::shared_ptr<ContPotential> _potential) :
    T(_T), potential(_potential) {};
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      boost::shared_ptr<Functor> funcExpEnergy (new ExpEnergy(T, potential) );
      boost::shared_ptr<Functor> funcPartition (new Partition(T, potential) );
      Maths::NumTech integrator;
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{	  
	  double energy;
	  if(i != steps.end() - 1)
	    {
	      it_step j = i + 1;
	      energy = 4.0 * M_PI / integrator.integrator(funcPartition, 
							  j->first, i->first, 100) *
		integrator.integrator(funcExpEnergy,j->first, i->first, 100);
	    }
	  else
	    {
	      energy = 4.0 * M_PI / integrator.integrator(funcPartition, 
							  1e-16, i->first, 100) *
		integrator.integrator(funcExpEnergy, 1e-16, i->first, 100);
	    }
	  i->second = energy;
	}
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
    double T;
    boost::shared_ptr<ContPotential> potential;
  };

  class Enr_Chapela: public Stepper
  {
  public:
    Enr_Chapela() {};
 
    virtual inline void genEnergy(std::vector<std::pair<double, double> >& steps)
    {
      if(steps.size() != 10)
	{
	  std::cerr << "Error: Chapela Energy stepping requires 10 steps" << std::endl;
	  exit(2);
	}
      double chapelaE[10] = {66.74, 27.55, 10.95, 3.81, 0.76, -0.47, -0.98, -0.55, -0.22, -0.06};
      int counter = 9;
      for(it_step i = steps.begin(); i != steps.end(); ++i)
	{
	  i->second = chapelaE[counter];
	  --counter;
	}
      
    }
  private:
    virtual inline void genPositions(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void genPotential(std::vector<std::pair<double, double> >& steps) {}
    virtual inline void addCore(std::vector<std::pair<double, double> >& steps) {}
  };
}
