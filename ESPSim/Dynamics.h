#pragma once

#include <vector>
#include <map>
#include <utility>
#include <iostream>

#include "Vector3.h"
#include "Particle.h"
#include "Stepmap.h"
#include "Simulator.h"

namespace Engine
{
  class Dynamics
  {
  public:
  Dynamics(Simulator* sim) :
    simulator(sim) {};
    void interact(double t, Scheduler::Event& event)
    {
      double mass = simulator->getParticles()[0].getMass();
      //update particle positions
      simulator->setParticles()[event.getP1()].move(t);
      simulator->setParticles()[event.getP2()].move(t);

      //calculate the separation properties
      PBCVector<double>r12(simulator->getSysLength(), true, 
			   simulator->getParticles()[event.getP1()].getR() 
			   - simulator->getParticles()[event.getP2()].getR());
      Vector3<double> v12 = simulator->getParticles()[event.getP1()].getV() 
	- simulator->getParticles()[event.getP2()].getV();
      double vdotr = v12.dotProd(r12.normalise());
      unsigned int step = simulator->getStepMap().getStep(event.getP1(), event.getP2());
      switch(event.getEventType())
	{
	case Scheduler::Event::IP_IN:
	  {
	    double dU = 0;
	    if(step == -1) //if no collision state found then particles must be outside outer step
	      dU = simulator->getSteps()[0].second; //energy is the outermost step height
	    else //if not outside then energy change is the difference in step heights
	      dU = simulator->getSteps()[step].second - simulator->getSteps()[step + 1].second;

	    if((vdotr * vdotr + 4.0 / mass * -dU) > 0) //if particles go over the step
	      {
		double A = -0.5 / mass * (vdotr + sqrt(vdotr * vdotr +  4.0 / mass * -dU)); //change in momentum
		//update particle velocities
		Vector3<double> deltav1 = (A / mass) * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		simulator->setStepMap().moveInwards(event.getP1(), event.getP2());
	      }
	    else //if bounce occurs
	      {
		Vector3<double> deltav1 = - vdotr * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
	      }
	    break;
	  }
	case Scheduler::Event::IP_OUT:
	  {
	
	    double dU = 0;
	    if(step == 0)
	      dU = -simulator->getSteps()[step].second;
	    else
	      dU = simulator->getSteps()[step - 1].second - simulator->getSteps()[step].second; 

	    if((vdotr * vdotr + 4.0 / mass * -dU) > 0) //if particles go over the step
	      {
		double A = -0.5 / mass * (vdotr - sqrt(vdotr * vdotr +  4.0 / mass * -dU)); //change in momentum

		//update particle velocities
		Vector3<double> deltav1 = A / mass * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
		simulator->setStepMap().moveOutwards(event.getP1(), event.getP2());
	      }
	    else //if bounce occurs
	      {
		Vector3<double> deltav1 = - vdotr * r12.normalise();
		simulator->setParticles()[event.getP1()].setV() += deltav1;
		simulator->setParticles()[event.getP2()].setV() -= deltav1;
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
