#include "Scheduler.h"
namespace Scheduler
{
  Event Scheduler::getNextEvent ()
  {
    size_t minIndex(0);
    double minTime(HUGE_VAL);
    for(size_t i(0); i < masterEL.size(); ++i)
      {
	if(masterEL[i].getCollisionTime() < minTime)
	  {
	    minTime = masterEL[i].getCollisionTime();
	    minIndex = i;
	  }
      }

    return masterEL[minIndex];    
  }
  void Scheduler::regenerate(double t, unsigned long long eventCount)
  {
    masterEL.clear();
    particleEL.clear(); 
    particleEL.resize(simulator->getParticles().size());
    sentPoint = simulator->getParticles().size();
    cellPoint = sentPoint + 1;
    BOOST_FOREACH(Particle p1, simulator->getParticles())
      {
	particleEL[p1.getID()].resize(cellPoint + 1);
	masterEL.push_back(getMinTime(t, p1.getID()));
      }
    
    if(masterEL.size() != simulator->getParticles().size())
      {
	std::cerr << "MasterEL not correct size" << std::endl;
	exit(3);
      }
    if(simulator->setThermostat()->is_initialised())
      {
	masterEL.push_back(Event());
	thermoPoint = masterEL.size() - 1;
	getThermoEvent(t, eventCount);
      }
  }
  void Scheduler::initialise()
  {
    regenerate(0,0);
  }

  void Scheduler::update(double t, unsigned int p1)
  {
    masterEL[p1] = getMinTime(t, p1);
  }
  void Scheduler::update(double t, unsigned int p1, unsigned int p2)
  {
    update(t, p1); update(t, p2);
  }

  Event Scheduler::getMinTime(double t, unsigned int p1)
  {
    nl->genNL(simulator->getParticles()[p1].getCell());
    for(std::vector<unsigned int>::iterator p2 = nl->getNeighbours().begin();
    	p2 != nl->getNeighbours().end(); ++p2)
      {
	if(p1 == *p2) continue; //if particle is itself move on
	Event::EventType eventType;
	double t_min_coll = getInteractionTime(p1, *p2, 
					       eventType);
	particleEL[p1][*p2] = Event(t + t_min_coll, p1, *p2, 
				    simulator->getParticles()[*p2].getNoColl(), eventType);
	//std::cerr << p1 << " & " << *p2 << " = " << t_min_coll + t << std::endl;
      }
    double sent_time = getSentinal(p1);
    particleEL[p1][sentPoint] = Event(t + sent_time, p1, -1, -1, Event::SENTINAL);
    double cell_time = nl->getCellTime(p1);
    particleEL[p1][sentPoint] = Event(t + cell_time, p1, -1, -1, Event::NEIGHBOURCELL);
    return *min_element(particleEL[p1].begin(), particleEL[p1].end());
  }

  double Scheduler::getSentinal(unsigned int p1)
  {
    double t_min = HUGE_VAL; //set the minimum time to infinity
    for(size_t dim (0); dim < 3; ++dim) //loop through all the dimensions
      {
	double vel = fabs(simulator->getParticles()[p1].getV()[dim]); //calculate the absolute velocity
	if(vel != 0)
	  {
	    double t_sent = 0.25 * (simulator->getSysLength() - simulator->getSteps()[0].first) / vel; //calculate the time when particle is invalid in that direction
	    t_min = (t_sent < t_min) ? t_sent : t_min; //if new calculated time is less than current minimum time
	  }
      }
    return t_min; 
  }
  
  void Scheduler::getThermoEvent(double t, unsigned long long eventCount)
    {
      
      Scheduler::masterEL[thermoPoint] 
	= Event(t + simulator->setThermostat()->getThermoTime(eventCount), 
		-1, -1, -1, Event::THERMOSTAT);
    }
  /*Event Scheduler::getRDF()
    {
    //sampler::createRDFEvent();
    //return Event(sampler::RDF::getTime(), -1, -1, Scheduler::Event::RDF);
    }
  */
  void Scheduler::invalidateEvent(double p1, double p2)
  {
    particleEL[p1][p2].setCollisionTime(HUGE_VAL);
  }
  double Scheduler::getInteractionTime(unsigned int p1, unsigned int p2, Event::EventType& eventType)
  {
    double t_min_out = HUGE_VAL, t_min_in = HUGE_VAL;
    //[[COMMENT]]can this be moved into a pair particles class...just to narrow this code down?
    PBCVector<double> r12(simulator->getSysLength(), 
			  true, 
			  simulator->getParticles()[p1].getR() 
			  - simulator->getParticles()[p2].getR());
    Vector3<double> v12 = simulator->getParticles()[p1].getV() 
      - simulator->getParticles()[p2].getV();
    //calculate dot products
    double r12sqr = r12.lengthSqr();
    double v12sqr = v12.lengthSqr();
    double vdotr = v12.dotProd(r12);
    int step = simulator->getStepMap().getStep(p1, p2);
    if (step == -1) //if particles are not in the step map then they must be outside range of influence, calculate when they come together
      {
	//calculate quadratic root arguments
	double c = r12sqr - simulator->getSteps()[0].first * simulator->getSteps()[0].first;
	double arg = vdotr * vdotr - v12sqr * c;
	if((vdotr < 0) && (arg >= 0)) //if particles come near enough to each other
	  t_min_in = c / (-vdotr + sqrt(arg)); //calculate time when particle's collide
	eventType = Event::IP_IN; //set the event type to an inward collision
	return (t_min_in > 0) ? t_min_in : HUGE_VAL; //if minimum time is positive then return time of collision 
      }
    else //if particles are in step map
      {
	if(step != simulator->getSteps().size() - 1) //if there is an innerstep to interact with
	  {
	    if(vdotr < 0) //if coming towards each other
	      {
		//calculate collision time inwards
		double c = r12sqr 
		  - simulator->getSteps()[step + 1].first * simulator->getSteps()[step + 1].first;
		double arg = vdotr * vdotr - v12sqr * c;
		if(arg >= 0) //if particles come near enough to each other
		  t_min_in = c / (-vdotr + sqrt(arg));
	      }
	  }

	//check for any outward steps
	double c = r12sqr - simulator->getSteps()[step].first * simulator->getSteps()[step].first;
	double arg = vdotr * vdotr - v12sqr * c;

	if(arg >= 0) //if particles come near enough to each other
	  t_min_out = (sqrt(arg) - vdotr) / v12sqr;
	else if(vdotr < 0) //There's been a numerical error, so just collide at their nearest point
	  t_min_out = -vdotr / v12sqr;

	if(t_min_out < -1)
	  {
	    std::cerr << "ERROR " << arg << " - " << vdotr << " - " << step
		      << std::endl; 
	    simulator->setStepMap().checkMap(simulator->getParticles(), simulator->getSteps(),
					     simulator->getSysLength());
	    exit(3);
	  }
	if(t_min_out < t_min_in || t_min_in < 0) //if outward time is smaller than inward time (or inward time is invalid)
	  {
	    //return event is an outward collision and return collision time
	    eventType = Event::IP_OUT;
	    return t_min_out;
	  }
	else //else use the inward collision time
	  {
	    eventType = Event::IP_IN;
	    return t_min_in;
	  }
      }
  }

}
