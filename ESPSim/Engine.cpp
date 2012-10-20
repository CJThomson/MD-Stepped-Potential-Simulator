#include "Engine.h"
namespace Engine
{
  void Engine::equilibrate()
  {
    std::ofstream equiLog;
    equiLog.open ("equiLog.dat"); //open file

    equilibration = true;
    std::cout << "\rEquilibration => Initialising Thermostat           " << std::flush;
    simulator->setThermostat()->initialise(simulator->getTemperature(), 
					   simulator->getRNG());
    std::cout << "\rEquilibration => Initialising Neighbout List       " << std::flush;
    boost::shared_ptr<NL::NL> nl = loadNL();
    //boost::shared_ptr<NL::NL> nl(new NL::NL_None());
    nl->initialise(simulator);
    std::cout << "\rEquilibration => Generating Event List             " << std::flush;
    Scheduler::Scheduler eventList(simulator, nl);
    eventList.initialise();
    Sampler::Sampler sampler(simulator->getParticles().size(),
			     simulator->getParticles()[0].getMass(),
			     simulator->getDensity(),false, false);
    sampler.initialise(simulator->getParticles(), simulator->getSteps(),
		       simulator->getStepMap());
    std::cout << "\rEquilibration => Starting Simulation               " << std::flush;
    while(simulator->isRunning(t, eventCount, true))
      {
	Scheduler::Event nextEvent= eventList.getNextEvent();
	handleEvent(nextEvent, eventList, sampler, nl);

	if(eventCount % 1000 == 0)
	  {
	    std::cout << "\rEquilibration => " << std::setprecision(4) 
		      << std::setprecision(3)<< std::setfill(' ') <<std::setw(5)
		      << simulator->progress(t, eventCount, true) * 100 << " %"
		      << std::setprecision(6) << std::setw(4)
		      << " t:" << t
		      << " U: " << sampler.getU() 
		      << " T: " << sampler.getT() << "      " 
		      << std::flush;
	  }
      }
    std::cout << "\rEquilibration => Complete                         " 
	      << "         " << std::flush;
    equiLog.close(); //close the file
  }
  void Engine::productionRun(bool firstRun, Logger::Logger& logger)
  {

    equilibration = false;
    if(firstRun)
      logger.init_Results(simulator->getSettings().getSampleColl());
    std::cout << "\rRunning => Initialising Thermostat                " << std::flush;
    simulator->setThermostat()->initialise(simulator->getTemperature(), 
					   simulator->getRNG());
    std::cout << "\rRunning => Initialising Neighbout List            " << std::flush;
    boost::shared_ptr<NL::NL> nl = loadNL();
    nl->initialise(simulator);
    std::cout << "\rRunning => Generating Event List                  " << std::flush;
    Scheduler::Scheduler eventList(simulator, nl);
    eventList.initialise();

    Sampler::Sampler sampler(simulator->getParticles().size(),
			     simulator->getParticles()[0].getMass(),
			     simulator->getDensity(), 
			     simulator->getSettings().getSampleColl(), false);
    sampler.initialise(simulator->getParticles(), simulator->getSteps(),
		       simulator->getStepMap());

    std::cout << "\rRunning => Starting Simulation                    " << std::flush;
    while(simulator->isRunning(t, eventCount, false))
      {
	Scheduler::Event nextEvent= eventList.getNextEvent();
	handleEvent(nextEvent, eventList, sampler, nl);
	if(eventCount % 1000 == 0)
	  {
	    if(eventCount > 2E6)
	      {
		std::cout << std::setfill(' ') <<std::setw(5) 
			  << "\rN:" << eventCount / 1e6 << "M";
	      }
	    else
	      {
		std::cout << std::setfill(' ') <<std::setw(5) 
			  << "\rN:" << eventCount / 1000 << "k";
	      }
	    std::cout << std::setfill(' ') << std::setw(5) 
		      << std::setprecision(3) << " t: " << t
		      << " T: " << sampler.getT() 
		      << " <T>: " << sampler.getMeanT()
		      << " U: " << sampler.getU() 
		      << " <U>: " << sampler.getMeanU() 
		      << " <P>: " << sampler.getP()
		      << std::flush;
	  }
      }
    std::cout << "\rRunning => Logger => Writing Results            " 
	      << "                  " << std::flush;
    sampler.terminate(eventCount, t);
    logger.update_Results(sampler,simulator->getSettings().getSampleColl());
    logger.write_Results(simulator->getSettings().getSampleColl());
    std::cout << "\rRunning => Complete                             " 
	      << "                  " << std::flush;
  }
  void Engine::handleEvent(Scheduler::Event& currentEvent, Scheduler::Scheduler& el,
			   Sampler::Sampler& sampler, boost::shared_ptr<NL::NL> nl)
  {
    switch (currentEvent.getEventType())
      {
      case Scheduler::Event::IP_IN:
      case Scheduler::Event::IP_OUT:
	{
	  handleInteraction(currentEvent, el, sampler);
	  break;
	}
      case Scheduler::Event::SENTINAL:
        {
	  handleSentinal(currentEvent, el, sampler);
	  break;
	}
      case Scheduler::Event::THERMOSTAT:
	{
	  freeStream(currentEvent.getCollisionTime(), sampler);
	  unsigned int particleID 
	    = simulator->setThermostat()->runThermostat(simulator->setParticles(), sampler);
	  simulator->setParticles()[particleID].incrNoColl();
	  el.update(t, particleID);
	  el.getThermoEvent(t, eventCount);
	  break;
	}
      case Scheduler::Event::NEIGHBOURCELL:
	{
	  freeStream(currentEvent.getCollisionTime(), sampler);
	  /*std::cerr << currentEvent.getP1() << " - " 
		    << currentEvent.getCollisionTime() << " - "
		    << simulator->getParticles()[currentEvent.getP1()].getCell() << " - "
		    << simulator->getParticles()[currentEvent.getP1()].getNextCell()
		    << std::endl;*/
	  nl->moveParticle(currentEvent.getP1());
	  el.update(t, currentEvent.getP1());
	  break;
	}
      case Scheduler::Event::RDF:
      case Scheduler::Event::NONE:
      default:
	std::cerr << "ERROR: Invalid event type encountered" << std::endl;
      }
  }

  void Engine::handleInteraction(Scheduler::Event& currentEvent, Scheduler::Scheduler& el,
				 Sampler::Sampler& sampler)
  {
    if(currentEvent.getP2Coll() != simulator->getParticles()[currentEvent.getP2()].getNoColl()) //check if collision is valid
      {
	/*std::cerr << "invalid collision" << std::endl; 
	std::cerr << "t = " << t << " pID = " << currentEvent.getP1() 
	  << " & " << currentEvent.getP2()<<std::endl;
	std::cerr << currentEvent.getP2Coll() << " - "
	<<  simulator->getParticles()[currentEvent.getP2()].getNoColl();*/
	el.update(t, currentEvent.getP1(), currentEvent.getP2());
      }
    else //if a valid event
      {
	freeStream(currentEvent.getCollisionTime(), sampler);
	/*std::cerr << "t = " << t << "pID = " << currentEvent.getP1() 
	  << " & " << currentEvent.getP2()<<std::endl;*/
	Dynamics dynamics(simulator);
	dynamics.interact(t, currentEvent, sampler);
	simulator->setParticles()[currentEvent.getP1()].incrNoColl();
	simulator->setParticles()[currentEvent.getP2()].incrNoColl();

	el.update(t, currentEvent.getP1(), currentEvent.getP2());
	++eventCount;
      }
  }
  void Engine::handleSentinal(Scheduler::Event& currentEvent, Scheduler::Scheduler& el,
			      Sampler::Sampler& sampler)
  {
    freeStream(currentEvent.getCollisionTime(), sampler);
    el.update(t, currentEvent.getP1());
  }
  void Engine::freeStream(double newT, Sampler::Sampler& sampler)
  {
    double dt = newT - t;
    for(std::vector<Particle>::iterator p = simulator->setParticles().begin();
	p != simulator->setParticles().end(); ++p)
      p->setR() += p->getV() * dt;
    if(!equilibration)
      sampler.freeStream(dt);
    t = newT; //update system time
  }
  boost::shared_ptr<NL::NL> Engine::loadNL()
  {
    boost::shared_ptr<NL::NL> nl(new NL::NL_None);
    if(strcmp(simulator->getSettings().getNLType(),"Simple") == 0)
      nl = boost::shared_ptr<NL::NL>(new NL::NL_Simple(simulator->getSteps()[0].first));
    return nl;
  }
} 
