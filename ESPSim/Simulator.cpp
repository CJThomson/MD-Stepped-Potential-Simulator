#include "Simulator.h"
#include "Engine.h"
void Simulator::loadSettings(int argc, char *argv[])
{
  //load the config file
  parseXML config(simSettings, simProperties, simPot, particles);

  namespace po = boost::program_options;
    
  //define all the program options classes
  po::options_description genericOpts("Generic options"), 
    simSetOpts("Simulator Settings"),
    simPropOpts("Simulation Properties"),
    simPotOpts("Potential Properties");

  //Add all the available options
  genericOpts.add_options()
    ("help", "produce help message")
    ("config,f", po::value<std::string>(),"config file to use")
    ;
    
  simSettings.setOptions(simSetOpts);
  simProperties.setOptions(simPropOpts);
  simPot.setOptions(simPotOpts);
  po::options_description cmdline_options;
  cmdline_options.add(genericOpts).add(simPropOpts).add(simSetOpts).add(simPotOpts);

  //put input into variables map
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  if(vm.count("help")) 
    {
      std::cout << cmdline_options << "\n";
      exit(0);
    }
  if(vm.count("config"))
      config.parseFile(vm["config"].as<std::string>().c_str());
  else
    config.parseFile(NULL);
  //process the variables map
  simSettings.loadCLSettings(vm);
  simProperties.loadCLSettings(vm);
  simPot.loadCLSettings(vm);
}

void Simulator::initialise()
{
  std::cout << "\rInitialisation => Random Number Generator      " << std::flush;
  RNG.seed();
  std::cout << "\rInitialisation => Particle => Positions        " << std::flush;

  particles.resize(simProperties.getN());
  boost::scoped_ptr<Lattice> initPos( new FCC);  
  initPos->placeParticles(particles, simProperties.getLength()); //initialise particle locations
  std::cout << "\rInitialisation => Particle => Velocity         " << std::flush;
  for(it_p p = particles.begin();p != particles.end(); ++p)
    p->resetV(RNG);

  zeroMomentum();
  std::cout << "\rInitialisation => Potential                    " << std::flush;
  initSteps();
  std::cout << "\rInitialisation => Pair Step Map                " << std::flush;
  stepmap.populateMap(particles, steps, simProperties.getLength());

  std::cout << "\rInitialisation => Pair Step Map => Checking    " << std::flush;
  stepmap.checkMap(particles, steps, simProperties.getLength());

  std::cout << "\rInitialisation => Logger => Create Out Config  " << std::flush;
  Logger::Logger logger;
  logger.write_outConfig(simSettings, simProperties, simPot);
  std::cout << "\rInitialisation => Complete                     " << std::endl;

}
void Simulator::equilibrate()
{
  std::cout << "\rEquilibration => Resetting Particles         " << std::flush;
  for(it_p p = particles.begin();p != particles.end(); ++p)
    p->reset();
  std::cout << "\rEquilibration => Initialising Engine         " << std::flush;
  Engine::Engine engine(this);
  std::cout << "\rEquilibration => Starting Simulation         " << std::flush;
  engine.equilibrate();
  std::cout << "\rEquilibration => Complete                    " << std::endl;

}

void Simulator::productionRun()
{
  Logger::Logger logger;
  for(size_t i(0); i < simSettings.getRuns(); ++i) //loop through several runs
    {
      std::cout << "\rRun Number: " << i+1 <<"                     " << std::endl;
      std::cout << "\rRunning => Resetting Particles               " << std::flush;
      for(it_p p = particles.begin();p != particles.end(); ++p)
	p->reset();
      
      std::cout << "\rRunning => Initialising Engine               " << std::flush;
      Engine::Engine engine(this);
      std::cout << "\rRunning => Starting Simulation               " << std::flush;
      engine.productionRun((i == 0 ? true : false ), logger);
      std::cout << "\rRunning => Complete                          " << std::flush;
    }
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
  if(simSettings.isTime(equilibration)) //if simulator is running for a certain length of time
    return (equilibration ? t < simSettings.getEQTime() : t < simSettings.getRunTime());
  else
    return (equilibration ? N < simSettings.getEQEvent() : N < simSettings.getRunEvent());
}

double Simulator::progress(double t, unsigned long long N, bool equilibration)
{
  if(simSettings.isTime(equilibration)) //if simulator is running for a certain length of time
    return (equilibration ? t / simSettings.getEQTime() : t / simSettings.getRunTime());
  else
    return (equilibration ? (double) N / simSettings.getEQEvent() 
	    : (double) N / simSettings.getRunEvent());
}
void Simulator::initSteps()
{
  boost::shared_ptr<Stepper::Stepper> stepper;
  bool genPotential = true;
  const char* option = simPot.getStepPotential();
  if(strcmp(option,"HardSphere") == 0)
    stepper.reset(new Stepper::Pot_HardCore(simPot.getSigma()));
  else if(strcmp(option, "SquareWell") == 0)
    stepper.reset(new Stepper::Pot_SquareWell(simPot.getSigma(), simPot.getEpsilon(),
					      simPot.getRCutOff()));
  else if(strcmp(option, "Chapela") == 0)
    stepper.reset(new Stepper::Pot_Chapela());
  else
    genPotential = false;

  if(genPotential)
    stepper->genPotential(steps);
  else
    {
      boost::shared_ptr<ContPotential> potential;
      option = simPot.getContPotential();
      if(strcmp(option,"LennardJones") == 0)
	potential.reset(new LennardJones(simPot.getSigma(), simPot.getEpsilon()));
      else if(strcmp(option,"LennardJones_shifted") == 0)
	potential.reset(new LennardJonesShifted(simPot.getSigma(), simPot.getEpsilon(),
						simPot.getRCutOff()));
      
      option = simPot.getStepPositions();
      if(strcmp(option, "Even") == 0)
	stepper.reset(new Stepper::Pos_Even(simPot.getRCutOff(), simPot.getNoStep()));
      else if(strcmp(option,"MeanForce") == 0)
	stepper.reset(new Stepper::Pos_ExpForce(simProperties.getT(),simPot.getRCutOff(), 
						simPot.getNoStep(), potential));
      else if(strcmp(option,"EvenEnergy") == 0)
	stepper.reset(new Stepper::Pos_EvenEnergy(simPot.getEnergyInterval(), 
						  simPot.getRCutOff(), potential));
      else if(strcmp(option,"Chapela") == 0)
	stepper.reset(new Stepper::Pos_Chapela());
      else 
	exit(0);
      stepper->genPositions(steps);
      option = simPot.getStepEnergies();
      if(strcmp(option, "Mid") == 0)
	stepper.reset(new Stepper::Enr_Mid(potential));
      else if (strcmp(option,"Virial") == 0 )
	stepper.reset(new Stepper::Enr_Virial(simProperties.getT(), potential));
      else if (strcmp(option,"Average") == 0)
	stepper.reset(new Stepper::Enr_Average(potential));
      else if (strcmp(option,"Chapela") == 0)
	stepper.reset(new Stepper::Enr_Chapela());
      stepper->genEnergy(steps);
    }
  

}
