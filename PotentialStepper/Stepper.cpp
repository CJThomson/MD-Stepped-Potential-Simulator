#include<math.h>
#include<iostream>
#include<vector>
#include<fstream>
#include "Stepper.h"
double Stepper::lj_eps = 1.0;
double Stepper::lj_sig = 1.0;
double Stepper::beta = 1.0 / 1.5;

int main()
{
  //double r_core = 0.966872;
  double rcutoff = 3.0;
  double noSteps = 10;
  std::vector<Steps> steps;
  Stepper stepper;

//  std::ofstream expectF;
//  expectF.open ("expectForce.dat"); //open file
//  for(double i = 0.001; i < 3.0; i+=0.001)
//    expectF << i << "\t" << stepper.generatePlot(i, 0.001) << std::endl; //output all diffusion coeffcients
//  expectF.close(); //close the file

  std::cerr << "finished" << std::endl;

  stepper.generateSteps(noSteps, rcutoff, Stepper::VIRIAL , Stepper::EXPECTEDFORCE, steps);
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
