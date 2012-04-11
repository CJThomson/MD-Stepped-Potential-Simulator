#include<math.h>
#include<iostream>
#include<vector>

double r_core = 0.1;
double r_cutoff = 3.0;
double noSteps = 10;
double beta = 1.0 / 1.34;
double lj_eps = 1.0;
double lj_sig = 1.0;
double partition_Function(double r);
double potential(double r);
double integrator_Simpson(double (*function)(double),
			  double upper_Bound, double lower_Bound, size_t bins);
double limit_solver(double target_area, double upper_bound, double lower_bound,
		    size_t integrator_intervals, size_t max_iterations, double tolerance);
int sign(double);

int main()
{
  double totalZ = integrator_Simpson(&partition_Function, r_cutoff, r_core, 1000);
  std::cout << totalZ << std::endl;
  std::vector<double> steps;
  double r_lower = r_core;
  for(size_t i(0); i < noSteps; ++i)
    {
      double step = limit_solver(totalZ / noSteps, r_cutoff, r_lower, 1e6, 100, 1e-7);
      std::cout << std::endl;
      if(step == 0)
	{
	  std::cerr << "ERROR: limit_solver could not find a step distance" << std::endl;
	  break;
	}
      else
	{
	  steps.push_back(step);
	  r_lower = step;
	}
    }
  for(int i = 0; i <steps.size(); ++i)
    {

      std::cout << "step: " << i  << " has width " << steps[i] << std::endl;
    }
} 
double partition_Function(double r)
{
  return 4 * M_PI * exp(-beta * potential(r)) * r * r;
}

double potential(double r)
{
  return  4.0 * lj_eps * (pow(lj_sig / r, 12) - pow(lj_sig / r, 6));
}
    
double integrator_Simpson(double (*function)(double),
			  double upper_Bound, double lower_Bound, size_t intervals)
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
	    
      sum += factor * function(lower_Bound + i * h);
    }
  return  h / 3.0 * sum;
}

double limit_solver(double target_area, double upper_bound, double lower_bound,
		    size_t integrator_intervals, size_t max_iterations, double tolerance)
{
  size_t iterations(0);
  double interval = (upper_bound - lower_bound) / 2.0;
  double b = lower_bound;

  while(iterations < max_iterations)
    {
      b += interval;
      double integral = integrator_Simpson(&partition_Function,
					   b, lower_bound, integrator_intervals);
      if(target_area - tolerance < integral &&  target_area + tolerance > integral)
	return b;

      else
       	std::cout << "\rlimit: " << b << " integral: " << integral << " target: " << target_area << std::flush;
      interval = fabs(interval / 2.0);
      if(integral > target_area)
	interval *= -1;
	
      ++iterations;
    }
    
  std::cerr << "ERROR: Maximum number of steps exceeded in limit_solver" << std::endl;
  return 0;
}

int sign(double x)
{
  return (x >= 0) ? 1 : -1;
}
