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
    VIRIAL,
    ENERGY,
    MID,
    AREA,
    ENERGYACTION
  } StepHeight;

  typedef enum {
    EVEN,
    PROBABILITY,
    EXPECTEDFORCE
  } StepWidth;
  //Constructors:
  // = = variables
  static double lj_eps; // lennard jones minimum energy
  static double lj_sig; // lennard jones distance of root
  static double beta; // inverse reduced temperature
  //Functions
  double generatePlot(double r, double r0)
  {
    return integrator_Simpson(&expected_Force, r, r0, 1000);
  }
  void generateSteps(unsigned int number_of_steps,
		     double r_cutoff,
		     StepHeight step_energy,
		     StepWidth step_radius,
		     std::vector<Steps>& genSteps)
    {
      genSteps.clear();

      //calculate the equivalent hard core
      //double r_core =  integrator_Simpson(&BHequivalentDiameter, lj_sig, ZERO, 1000);
      double r_core = ZERO;
      //genSteps.push_back(Steps(r_core,0));
      totalZ = integrator_Simpson(&partition_Function, r_cutoff, r_core, 1000);
      //--number_of_steps;
      switch(step_radius)
	{
	case EVEN:
	  {
	    for(size_t i(1); i <= number_of_steps; ++i) //generate step lengths
	      {
		double step = i * (r_cutoff - r_core) / number_of_steps + r_core;
		genSteps.push_back(Steps(step,0));
	      }
	    break;
	  }
	case PROBABILITY:
	  {
	    //calculate total partition funciton

	    double r_lower = 0.9;
	    for(size_t i(0); i < number_of_steps; ++i) //generate step lengths
	      {
		double step = limit_solver(&partition_Function, totalZ / number_of_steps, r_cutoff, r_lower, 100, 1e6, 1e-10);
		if(step == 0)
		  break;
		else
		  {
		    genSteps.push_back(Steps(step,0));
		    r_lower = step;
		  }
	      }
	  }
	  break;
	case EXPECTEDFORCE:
	  //calculate total partition funciton
	  double totalEF = integrator_Simpson(&expected_Force,r_cutoff,  r_core, 1000);
	  double r_lower = r_core;
	  for(size_t i(0); i < number_of_steps - 1; ++i) //generate step lengths
	    {
	      double step = limit_solver_bisection(&expected_Force, totalEF / number_of_steps, r_cutoff, r_lower, 1000, 1e6, 1e-5);
	      if(step == 0)
		break;
	      else
		{
		  genSteps.push_back(Steps(step,0));
		  r_lower = step;
		}
	    }
	  genSteps.push_back(Steps(r_cutoff,0));
	  break;
	}


      for(std::vector<Steps>::iterator step_i = genSteps.begin();
	  step_i != genSteps.end(); ++step_i) //generate step lengths
	{
	  double energy = 0;
	  switch(step_energy)
	    {
	    case VIRIAL:
	      if(step_i->step_energy == 0)
		{
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
		}
	      break;
	    case ENERGY:
	      if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  energy = 4.0 * M_PI
		    * 1.0 / integrator_Simpson(&partition_Function,
					       step_i->step_radius,
					       step_j->step_radius,
					       100)
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
					       100)
		    * integrator_Simpson(&internal_Energy,
					 step_i->step_radius,
					 ZERO,
					 100);
		}
	      step_i->step_energy = energy;
	      break;
	    case MID:

	      if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  energy = potential((step_i->step_radius + step_j->step_radius) * 0.5);
		}
	      else
		{
		  energy = potential(step_i->step_radius);
		}
	      step_i->step_energy = energy;
	      break;
	    case AREA:
	       if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  double area = integrator_Simpson(&potential, step_i->step_radius, step_j->step_radius, 100);
		  energy = area / (step_i->step_radius - step_j->step_radius);
		}
	      else
		{
		  energy = 40;
		}
	      step_i->step_energy = energy;
	      break;
	    case ENERGYACTION:
	      if(step_i != genSteps.begin())
		{
		  std::vector<Steps>::iterator step_j = step_i - 1;
		  energy = expected_Force((step_i->step_radius + step_j->step_radius) * 0.5);
		}
	      else
		{
		  energy = expected_Force(step_i->step_radius);
		}
	      step_i->step_energy = energy;
	      break;
	    }
	}
    }


 private:
  double totalZ;
  // = = Functions
  // = Function to calcualte Baker Henderson Equivalent Hard Sphere Diameter =
  static inline double BHequivalentDiameter(double z)
  {
    return (1.0 - exp(-beta * potential(z)));
  }
  // = Expected force function thing
  static inline double expected_Force(double r)
  {
    return fabs(force(r) * partition_Function(r));
  }
  // = Partition Function
  static inline double partition_Function(double r)
  {
    return 4 * M_PI * exp(-beta * potential(r)) * r * r;
  }
  // = Internal Energy Integral
  static inline double internal_Energy(double r)
  {
    return potential(r) * exp(-beta * potential(r)) * r * r;
  }
  // = Lennard Jones Potential
  static inline double potential(double r)
  {
    return 4.0 * lj_eps * (pow(lj_sig / r, 12) - pow(lj_sig / r, 6));
  }
  // = Lennard Jones Force
  static inline double force(double r)
  {
    return 24.0 * lj_eps / lj_sig * (2.0 * pow(lj_sig / r, 13) - pow(lj_sig / r, 7));
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
  double limit_solver(double (*function)(double), double target_area,
		      double upper_bound, double lower_bound,
		      size_t integrator_intervals, size_t max_iterations, double tolerance)
  {
    size_t iterations(0);
    double interval = (upper_bound - lower_bound) / 2.0;
    double b = upper_bound;
    while(iterations < max_iterations)
      {

	double integral = integrator_Simpson(function,
						  b, lower_bound, integrator_intervals);
	double root = integral - target_area;
	std::cerr << b << " - " << integral << " - " << target_area << std::endl; 
	if(root - tolerance < 0 &&  root + tolerance > 0)
	  return b;

	double differential = function(b);
	b -= root / differential;
	++iterations;
      }
    std::cerr << "ERROR: Maximum number of steps exceeded in limit_solver" << std::endl;
    return 0;
    }
  // = Function to calculate limits of integrals
  double limit_solver_bisection(double (*function)(double), double target_area,
				double upper_bound, double lower_bound,
				size_t integrator_intervals, size_t max_iterations, double tolerance)
  {
    size_t iterations(0);
    double interval = upper_bound - lower_bound;
    double b = lower_bound;
    while(iterations < max_iterations)
      {
	interval *= 0.5;
	double integral = integrator_Simpson(function,
					     b + interval, lower_bound, integrator_intervals);
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
