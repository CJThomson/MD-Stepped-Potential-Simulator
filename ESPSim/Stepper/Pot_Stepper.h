#pragma once
#include <vector>
#include <utility>
#include <limits>

#include <boost/scoped_ptr.hpp>
#include "Stepper.h"
#include "Pos_Stepper.h"
#include "Enr_Stepper.h"

namespace Stepper
{
  class Pot_HardCore: public Stepper
  {
  public:
  Pot_HardCore(double _sigma) :
    sigma(_sigma) {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      steps.push_back(std::make_pair(sigma,std::numeric_limits<double>::max()));
    }
  private:
    double sigma;
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };

  class Pot_SquareWell: public Stepper
  {
  public:
  Pot_SquareWell(double _sigma, double _epsilon, double _lambda) : 
    sigma(_sigma), epsilon(_epsilon), lambda(_lambda) {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      steps.push_back(std::make_pair(lambda * sigma, epsilon)); //the well
      steps.push_back(std::make_pair(sigma, std::numeric_limits<double>::max())); //the core
    }
  private:
    double sigma;
    double epsilon;
    double lambda;
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };


  class Pot_Chapela: public Stepper
  {
  public:
    Pot_Chapela() {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      {
	boost::scoped_ptr<Stepper> chapelaPos(new Pos_Chapela());
	chapelaPos->genPositions(steps);
      }
      {
	boost::scoped_ptr<Stepper> chapelaEnr(new Enr_Chapela());
	chapelaEnr->genEnergy(steps);
      }
    }
  private:
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };
}
