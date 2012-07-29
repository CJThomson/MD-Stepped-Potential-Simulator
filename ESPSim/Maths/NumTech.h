#pragma once
#include <iostream>

#include "Functions.h"
namespace Maths
{
  class NumTech
  {

  public:
    double integrator(Functions* func,
		      double lower_Bound, double upper_Bound, int intervals)
    {
      if(intervals % 2 != 0)
	{
	  std::cerr << "ERROR: Number of intervals for Simpson's Rule must be even" << std::endl;
	  return 0;
	}
      double h = (upper_Bound - lower_Bound) / intervals;
      double sum = 0;
      for (size_t i (0); i <= intervals; ++i)
	{
	  double factor = 1.0;
	  if(i > 0 && i < intervals)
	    factor = (i % 2 == 0) ? 2.0 : 4.0;
	  sum += factor * func->f(lower_Bound + i * h);
	}
      return  h / 3.0 * sum;
    }

    // = Function to calculate limits of integrals
    double limitSolver(Functions* func, double target_area,
		       double lower_bound, double upper_bound,
		       size_t integrator_intervals, size_t max_iterations, double tolerance)
    {
      size_t iterations(0);
      double interval = upper_bound - lower_bound;
      double b = lower_bound;
      while(iterations < max_iterations)
	{
	  interval *= 0.5;
	  double integral = integrator(func, b + interval, lower_bound, integrator_intervals);
	  double root = integral - target_area;
	  if(root - tolerance < 0 &&  root + tolerance > 0)
	    return b + interval;

	  if(integral < target_area)
	    b += interval;
	  ++iterations;
	}
      std::cerr << "ERROR: Maximum number of steps exceeded in limit_solver" << std::endl;
      return 0;
    }



  };
}
