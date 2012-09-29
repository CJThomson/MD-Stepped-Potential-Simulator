#include "Engine.h"
namespace Engine
{
  void Engine::equilibrate()
  {
    std::ofstream equiLog;
    equiLog.open ("equiLog.dat"); //open file

    equilibration = true;
    std::cout << "\rEquilibration => Generating Event List             " << std::flush;
    Scheduler::Scheduler eventList(simulator);
    eventList.initialise();
    Sampler::Sampler sampler(simulator->getParticles().size(),
			     simulator->getParticles()[0].getMass(),
			     simulator->getDensity(),false, false);
    sampler.initialise(simulator->getParticles(), simulator->getSteps(),
		       simulator->getStepMap());
    std::cout << "\rEquilibration => Starting Simulation               " << std::flush;
    while(simulator->isRunning(t, eventCount, true))
      {
	/*simulator->setStepMap().checkMap(simulator->getParticles (), simulator->getSteps(),
					 simulator->getSysLength());
					 //eventList.initialise();*/
	Scheduler::Event nextEvent= eventList.getNextEvent();
	handleEvent(nextEvent, eventList, sampler);

	if(int(simulator->progress(t, eventCount, true) * 10000) % 50 == 0)
	  {
	    std::cout << "\rEquilibration => " << std::setprecision(4) 
		      << std::setprecision(3)<< std::setfill(' ') <<std::setw(5)
		      << simulator->progress(t, eventCount, true) * 100 << " %"
		      << std::setprecision(6) << std::setw(4)
		      << " K:" << sampler.getKE() / simulator->getParticles().size()
		      << " U: " << sampler.getU() 
		      << " T: " << sampler.getT() << "      " 
		      << std::flush;
	  }
      }
    std::cout << "\rEquilibration => Complete                         " 
	      << "         " << std::flush;
    equiLog.close(); //close the file
  }
  void Engine::productionRun()
  {
    equilibration = false;
    std::cout << "\rRunning => Generating Event List                  " << std::flush;
    Scheduler::Scheduler eventList(simulator);
    eventList.initialise();

    Sampler::Sampler sampler(simulator->getParticles().size(),
			     simulator->getParticles()[0].getMass(),
			     simulator->getDensity(),false, false);
    sampler.initialise(simulator->getParticles(), simulator->getSteps(),
		       simulator->getStepMap());

    std::cout << "\rRunning => Starting Simulation                    " << std::flush;
    while(simulator->isRunning(t, eventCount, false))
      {
	Scheduler::Event nextEvent= eventList.getNextEvent();
	handleEvent(nextEvent, eventList, sampler);
	if(!simulator->setStepMap().checkMap(simulator->getParticles (), simulator->getSteps(),
					     simulator->getSysLength()))
	  eventList.initialise();
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
		      << " <P>: " << sampler.getP() << "-" << sampler.getMomFlux()
		      << "-" << sampler.getFluxTime()
		      << std::flush;
	  }
      }
    std::cout << "\rEquilibration => Complete                    " 
	      << "         " << std::flush;
  }
  void Engine::handleEvent(Scheduler::Event& currentEvent, Scheduler::Scheduler& el,
			   Sampler::Sampler& sampler)
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
      case Scheduler::Event::NEIGHBOURCELL:
      case Scheduler::Event::THERMOSTAT:
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
	el.update(t, currentEvent.getP1(), currentEvent.getP2());
    else //if a valid event
      {
	freeStream(currentEvent.getCollisionTime(), sampler);
	Dynamics dynamics(simulator);
	dynamics.interact(t, currentEvent, sampler);
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
    if(!equilibration)
      sampler.freeStream(dt);
    t = newT; //update system time
  }
} 
