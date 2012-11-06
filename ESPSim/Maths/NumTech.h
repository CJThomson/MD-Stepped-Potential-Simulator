#pragma once
#include <iostream>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include "Functor.h"
namespace Maths
{
  //a template sign function
  template < typename T> int sgn(T x)
    {
      if(x < 0) return -1; 
      if(x > 0) return 1; 
      return 0;
    }

  class NumTech
  {

  public:
    double integrator(boost::shared_ptr<Functor> f,
		      double lower_Bound, double upper_Bound, 
		      unsigned int intervals)
    {
      //Carries out an integration using: = = Simpson's Rule = =
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
	  sum += factor * (*f)(lower_Bound + i * h);
	}
      return  h / 3.0 * sum;
    }

    double rootFinder(boost::shared_ptr<Functor> f, 
		      double lower_bound, double upper_bound,
		      unsigned int max_iterations, double tolerance,
		      double target = 0)
    {
      //Solves for the roots of an equation using: = = Bisection = =
      double midpoint, f_mid;
      double f_lower = (*f)(lower_bound);
      for(size_t iterations = 0; iterations < max_iterations; ++iterations)
	{
	  //calculate the mid point
	  midpoint = 0.5 * ( upper_bound + lower_bound );
	  f_mid = (*f)(midpoint) - target;
	  if( abs(f_mid) < tolerance ) //if error is acceptable
	    return midpoint;
	  //if the sign of mid and lower are the same, root must be
	  //between mid and upper
	  if(sgn(f_mid) == sgn(f_lower))
	    {
	      lower_bound = midpoint;
	      f_lower = f_mid;
	    }
	  else
	    upper_bound = midpoint;
	}
      std::cerr << "ERROR: Maximum number of steps exceeded in limit_solver" << std::endl;
      return 0;
    }

  };
}
