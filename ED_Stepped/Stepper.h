#pragma once
/* Class to produce steps model a soft potential */
struct Steps
{
Steps(double r, double h) : step_radius(r), step_energy(h) {}
  double step_energy;
  double step_radius;
};

class Stepper
{
#define ZERO 1e-16
 public:
  typedef enum {
    CHAPELA,
    PROBABILITY,
    ENERGY
  } StepType;

  //Constructors:
  // = = variables
  static double lj_eps; // lennard jones minimum energy
  static double lj_sig; // lennard jones distance of root
  static double beta; // inverse reduced temperature
  //Functions
  void generateSteps(unsigned int number_of_steps,
		     double r_cutoff,
		     StepType step_type,
		     std::vector<Steps>& genSteps)
    {
      genSteps.clear();

      //calculate the equivalent hard core
      //double r_core =  integrator_Simpson(&BHequivalentDiameter, lj_sig, ZERO, 1000);
      double r_core = 0.8;
      genSteps.push_back(Steps(r_core, 0));
      //calculate total partition funciton
      double totalZ = integrator_Simpson(&partition_Function, r_cutoff, r_core, 1000);
      double r_lower = r_core;
      for(size_t i(0); i < number_of_steps; ++i) //generate step lengths
	{
	  double step = limit_solver(totalZ / number_of_steps, r_cutoff, r_lower, 100, 1e6, 1e-10);
	  if(step == 0)
	    {
	      break;
	    }
	  else
	    {
	      genSteps.push_back(Steps(step,0));
	      r_lower = step;
	    }
	}

      for(std::vector<Steps>::iterator step_i = genSteps.begin();
	  step_i != genSteps.end(); ++step_i) //generate step lengths
	{
	  double energy = 0;
	  switch(step_type)
	    {
	    case PROBABILITY:

	      if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  energy = - 1.0 / beta
		    * log(3.0 / ( 4.0 * M_PI * (pow(step_i->step_radius, 3) - pow((*(step_i - 1)).step_radius, 3)))
			  * integrator_Simpson(&partition_Function,
					       step_i->step_radius,
					       step_j->step_radius,
					       100));
		}
	      else
		{
		  energy = - 1.0 / beta
		    * log(3.0 / ( 4.0 * M_PI * (pow(step_i->step_radius, 3)))
			  * integrator_Simpson(&partition_Function,
					     step_i->step_radius,
					     ZERO,
					     100));
		}
	      step_i->step_energy = energy;
	      break;
	    case ENERGY:

	      if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  energy = 4.0 * M_PI
		    * 1.0 / integrator_Simpson(&partition_Function,
					       step_i->step_radius,
					       step_j->step_radius,
					       1000)
		    * integrator_Simpson(&internal_Energy,
					 step_i->step_radius,
					 step_j->step_radius,
					 100);
		}
	      else
		{
		  energy = 4.0 * M_PI
		    * 1.0 / integrator_Simpson(&partition_Function,
					       step_i->step_radius,
					       ZERO,
					       1000)
		    * integrator_Simpson(&internal_Energy,
					 step_i->step_radius,
					 ZERO,
					 100);
		}
	      step_i->step_energy = energy;
	      break;
	    }
	}
    }
 private:

  // = = Functions
  // = Function to calcualte Baker Henderson Equivalent Hard Sphere Diameter =
   static inline double BHequivalentDiameter(double z)
   {
     return (1.0 - exp(-beta * potential(z)));
   }
  // = Partition Function
   static inline double partition_Function(double r)
   {
     return 4 * M_PI * (1 - exp(-beta * potential(r))) * r * r;
   }
  // = Internal Energy Integral
   static inline double internal_Energy(double r)
   {
    return potential(r) * exp(-beta * potential(r)) * r * r;
   }
  // = Lennard Jones Potential
   static inline double potential(double r)
   {
     return  4.0 * lj_eps * (pow(lj_sig / r, 12) - pow(lj_sig / r, 6));
   }

  // = Numerical Integrator using Simpson's Rule
   double integrator_Simpson(double (*function)(double),
			    double upper_Bound, double lower_Bound, int intervals)
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
	Stepper obj;
	sum += factor * function(lower_Bound + i * h);
      }
    return  h / 3.0 * sum;
  }

  // = Function to calculate limits of integrals
  double limit_solver(double target_area, double upper_bound, double lower_bound,
		      size_t integrator_intervals, size_t max_iterations, double tolerance)
  {
    size_t iterations(0);
    double interval = (upper_bound - lower_bound) / 2.0;
    double b = upper_bound;
    while(iterations < max_iterations)
      {
	double integral = integrator_Simpson(&partition_Function,
						  b, lower_bound, integrator_intervals);
	double root = integral - target_area;
	if(root - tolerance < 0 &&  root + tolerance > 0)
	  return b;

	double differential = partition_Function(b);
	b -= root / differential;
	++iterations;
      }
    std::cerr << "ERROR: Maximum number of steps exceeded in limit_solver" << std::endl;
    return 0;
    }
};
