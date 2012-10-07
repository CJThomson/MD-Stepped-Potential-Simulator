#pragma once

#include <vector>
#include <map>
#include <utility>
#include <iostream>

#include "Vector3.h"
#include "Particle.h"
#include "Stepmap.h"
#include "Sampler.h"

namespace Engine
{
  class Dynamics
  {
  public:
  Dynamics(Simulator* sim) :
    simulator(sim) {};
    void interact(double t, Scheduler::Event& event, Sampler::Sampler& sampler)
    {
      //get the masses of the particles
      double p1mass = simulator->getParticles()[event.getP1()].getMass();
      double p2mass = simulator->getParticles()[event.getP2()].getMass();
      double mass = 2 * p1mass * p2mass / (p1mass + p2mass);//fix to correct mu

      double dU = 0;

      //update particle positions
      //calculate the separation properties
      PBCVector<double>r12(simulator->getSysLength(), true, 
			   simulator->getParticles()[event.getP1()].getR() 
			   - simulator->getParticles()[event.getP2()].getR());
      Vector3<double> v12 = simulator->getParticles()[event.getP1()].getV() 
	- simulator->getParticles()[event.getP2()].getV();
      double vdotr = v12.dotProd(r12.normalise());
      Stepmap::it_map it_step = simulator->setStepMap().setStepPntr(event.getP1(), event.getP2());
      switch(event.getEventType())
	{
	case Scheduler::Event::IP_IN:
	  {
	    if(it_step == simulator->setStepMap().getEndPntr()) //if no collision state found then particles must be outside outer step
	      dU = simulator->getSteps()[0].second; //energy is the outermost step height
	    else //if not outside then energy change is the difference in step heights
	      dU = simulator->getSteps()[it_step->second + 1].second 
		- simulator->getSteps()[it_step->second].second;

	    if((vdotr * vdotr - 4.0 / mass * dU) > 0) //if particles go over the step
	      {
		double A = -0.5 / mass * (vdotr + sqrt(vdotr * vdotr - 4.0 / mass * dU)); //change in momentum
		//update particle velocities
		Vector3<double> deltav1 = (A / mass) * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		sampler.changeMomentumFlux(r12.dotProd(deltav1));
		sampler.changePotential(dU);
		sampler.eventCount(2, (it_step == simulator->setStepMap().getEndPntr() 
				       ? -1 : it_step->second), true);
		simulator->setStepMap().moveInwards(event.getP1(), event.getP2());
	      }
	    else //if bounce occurs
	      {
		Vector3<double> deltav1 = - vdotr * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		sampler.changeMomentumFlux(r12.dotProd(deltav1));
		sampler.eventCount(4, (it_step == simulator->setStepMap().getEndPntr() 
				       ? -1 : it_step->second), true);
	      }
	    break;
	  }
	case Scheduler::Event::IP_OUT:
	  {

	    if(it_step->second == 0)
	      dU = -simulator->getSteps()[it_step->second].second;
	    else
	      dU = simulator->getSteps()[it_step->second - 1].second 
		- simulator->getSteps()[it_step->second].second; 

	    if((vdotr * vdotr - 4.0 / mass * dU) > 0) //if particles go over the step
	      {
		double A = -0.5 / mass * (vdotr - sqrt(vdotr * vdotr - 4.0 / mass * dU)); //change in momentum

		//update particle velocities
		Vector3<double> deltav1 = A / mass * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		sampler.changeMomentumFlux(r12.dotProd(deltav1));
		sampler.changePotential(dU);
		sampler.eventCount(2, it_step->second, false);
		simulator->setStepMap().moveOutwards(event.getP1(), event.getP2()); 
	      }
	    else //if bounce occurs
	      {
		Vector3<double> deltav1 = - vdotr * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		sampler.changeMomentumFlux(r12.dotProd(deltav1));
		sampler.eventCount(4, it_step->second, false);
	      }
	    break;
	  }
	default:
	  std::cerr << "Unhandled calcVelocities type" << std::endl;
	  std::exit(1);
	}

    }

  private:
    Simulator* simulator;
  };
  
}
