#pragma once
#include <vector>
#include <utility>
#include <limits>

#include "Stepper.h"

namespace Stepper
{
  class Pot_HardCore: public Stepper
  {
  public:
  Pot_HardCore(double sigma) :
    Stepper(sigma, 0, 0, 0, 0) {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      steps.push_back(std::make_pair(sigma,std::numeric_limits<double>::max()));
    }
  private:
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };

  class Pot_SquareWell: public Stepper
  {
  public:
  Pot_SquareWell(double sigma, double epsilon, double lambda) : 
    Stepper(sigma, epsilon, lambda, 0, 0) {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      steps.push_back(std::make_pair(lambda * sigma, epsilon)); //the well
      steps.push_back(std::make_pair(sigma, std::numeric_limits<double>::max())); //the core
    }
  private:
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };


  class Pot_Chapela: public Stepper
  {
  public:
  Pot_Chapela() :
    Stepper(0, 0, 0, 0, 0) {};
    virtual void genPotential(std::vector<std::pair<double, double> >& steps)
    {
      steps.clear();
      double chapelaR[10] = {0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.25, 1.45, 1.75, 2.3};
      double chapelaE[10] = {66.74, 27.55, 10.95, 3.81, 0.76, -0.47, -0.98, -0.55, -0.22, -0.06};
      for(size_t i = 0; i < 10; ++i)
	steps.push_back(std::make_pair(chapelaR[9 - i], chapelaE[9 - i]));
    }
  private:
    virtual void addCore(std::vector<std::pair<double, double> >& steps) {};
    virtual void genPositions(std::vector<std::pair<double, double> >& steps) {};
    virtual void genEnergy(std::vector<std::pair<double, double> >& steps) {};
  };
}
