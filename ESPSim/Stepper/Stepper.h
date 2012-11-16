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

  };

}

