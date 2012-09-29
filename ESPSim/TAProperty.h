#pragma once

namespace Sampler
{
  class TAProperty
  {
  public:
    TAProperty() { curValue = t = average = stdDev = 0; }
    TAProperty(const double& value) { TAProperty(); curValue = value; }
    TAProperty& operator= (const double& value) { curValue = value; return *this; }
    TAProperty& operator+= (const double& value) { curValue += value; return *this; }
    TAProperty& operator-= (const double& value) { curValue -= value; return *this; }
    void stream (const double dt) 
    { 
      t += dt;
      average += current() * dt;
      stdDev += current() * current() * dt;
    }

    double current() const { return curValue; }
    double mean() const { return t == 0 ? curValue : average / t; }
    double meanSqr() const { return t == 0 ? curValue * curValue : stdDev / t; }
    double time() const { return t; }
  private:
    double t;
    double curValue;
    double average;
    double stdDev;
  };

}
