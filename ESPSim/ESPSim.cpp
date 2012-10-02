/* This is my new stepped potential simualator born out of the burned
out remains of my previous Soft_Event sim */

//These should go in a declares header eventually
#include <iostream>
#include <boost/program_options.hpp>

#include "Simulator.h"
using namespace std;

int main(int argc, char *argv[])
{
  cout << "=  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =" << endl;
  cout << "=                                                           =" << endl;
  cout << "=   ESPSim - The Event-driven Stepped Potential Simulator   =" << endl;
  cout << "=                   Chris J Thomson (2012)                  =" << endl;
  cout << "=                                                           =" << endl;
  cout << "=  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =" << endl;
  cout << endl;
  Simulator sim;

  sim.loadSettings(argc, argv);  //load all the system settings
  sim.initialise(); //initialise the system
  sim.equilibrate(); //equilibrate the system
  sim.productionRun(); //run the production runs
  cout << endl;
  return 0;
}

