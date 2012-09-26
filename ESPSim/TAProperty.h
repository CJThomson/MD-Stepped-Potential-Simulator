#pragma once

namespace Sampler
{
  class TAProperty
  {
  public:
    TAProperty() { curValue = 0; }
    TAProperty(const double& value) { curValue = value; }
    TAProperty& operator= (const double& value) { curValue = value; return *this; }
    TAProperty& operator+= (const double& value) { operator=(curValue + value); }
    TAProperty& operator-= (const double& value) { operator=(curValue - value); }
    void stream (const double& dt) 
    { 
      t += dt;
      average += current() * dt;
      stdDev += current() * current() * dt;
    }

    double current() const { return curValue; }
    double mean() const { return t == 0 ? current() : average / t; }
    double meanSqr() const { return t == 0 ? current() * current() : stdDev / t; }
    double time() const { return t; }
  private:
    double t;
    double curValue;
    double average;
    double stdDev;
  };

}
