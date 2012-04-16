#include<math.h>
#include<iostream>
#include<vector>

#include "Stepper.h"
double Stepper::lj_eps = 1.0;
double Stepper::lj_sig = 1.0;
double Stepper::beta = 1.0 / 1.34;

int main()
{
  //double r_core = 0.966872;
  double rcutoff = 2.3;
  double noSteps = 9;


  std::vector<Steps> steps;
  Stepper stepper;

  stepper.generateSteps(noSteps, rcutoff, Stepper::PROBABILITY, steps);
  int counter = 0;
  for(std::vector<Steps>::iterator step_i = steps.begin();
      step_i != steps.end(); ++step_i) //generate step lengths
    {
      
      std::cout << "Step " << counter 
		<< " has radius " << step_i->step_radius 
		<< " and energy " << step_i->step_energy << std::endl;
      ++counter;
    }
}
