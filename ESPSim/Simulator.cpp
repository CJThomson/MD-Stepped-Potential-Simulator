#include "Simulator.h"

void Simulator::loadSettings(int argc, char *argv[])
{
    
  //Then code the config file/cmdline using boost
  namespace po = boost::program_options;
    
  //define all the program options classes
  po::options_description genericOpts("Generic options"), 
    simSetOpts("Simulator Settings"),
    simPropOpts("Simulation Properties");

  //Add all the available options
  genericOpts.add_options()
    ("help", "produce help message")
    ;
    
  simSettings.setOptions(simSetOpts);
  simProperties.setOptions(simPropOpts);
    
  po::options_description cmdline_options;
  cmdline_options.add(genericOpts).add(simPropOpts).add(simSetOpts);

  po::options_description config_file_options;
  config_file_options.add(simSetOpts).add(simPropOpts);

  //put input into variables map
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);

  //process the variables map
  simSettings.loadCLSettings(vm);
  simProperties.loadCLSettings(vm);
}

void Simulator::initialise()
{
  std::cout << "\rInitialising => Random Number Generator    " << std::flush;
  RNG.seed();
  std::cout << "\rInitialising => Particle => Positions      " << std::flush;

  particles.resize(simProperties.getN());
  Lattice* initPos = new FCC; 
  {
    initPos->placeParticles(particles, simProperties.getLength()); //initialise particle locations
    delete initPos; //release memory
  }
  std::cout << "\rInitialising => Particle => Velocity       " << std::flush;
  for(it_p p = particles.begin();p != particles.end(); ++p)
    p->resetV(RNG);
  zeroMomentum();
  std::cout << "\rInitialising => Potential                  " << std::flush;
  ContPotential* potential = new LennardJones(1,1);
  {
    Stepper::Stepper* stepper = new Stepper::Pos_Even(3.0, 10);
    {
      stepper->genPositions(steps);
      delete stepper;
    }
    stepper = new Stepper::Enr_Mid(potential);
    {
      stepper->genEnergy(steps);
      delete stepper;
    }
    delete potential;
  }
  std::cout << "\rInitialising => Pair Step Map              " << std::flush;
  stepmap.populateMap(particles, steps, simProperties.getLength());
  std::cout << "\rInitialising => Pair Step Map => Checking  " << std::flush;
  stepmap.checkMap(particles, steps, simProperties.getLength());
  std::cout << "\rInitialising => Complete                   " << std::flush;

}
void Simulator::equilbriate()
{
  std::cout << "\rEquilibration => Resetting Particles       " << std::flush;
  for(it_p p = particles.begin();p != particles.end(); ++p)
    p->reset();
  std::cout << "\rEquilibration => Initialising Engine       " << std::flush;
  Engine engine(this, &particles);
  std::cout << "\rEquilibration => Starting Simulation       " << std::flush;
    
  std::cout << "\rEquilibration => Complete                  " << std::flush;

}
void Simulator::zeroMomentum()
{
  Vector3<double> sum;
  for(it_p p = particles.begin(); p != particles.end(); ++p) 
    sum += p->getV();
  
  for(it_p p = particles.begin(); p != particles.end(); ++p) 
    p->setV() -= sum / particles.size();
}

bool Simulator::isRunning(double t, unsigned long long N, bool equilibration)
{
  if(simSettings.isTime()) //if simulator is running for a certain length of time
    return (equilibration ? t < simSettings.getEQTime() : t < simSettings.getRunTime());
  else
    return (equilibration ? N < simSettings.getEQEvent() : N < simSettings.getRunEvent());
}
