#pragma once

#include <vector>
#include <iomanip>
#include <iostream>

#include "Particle.h"
#include "Simulator.h"
#include "Scheduler.h"
#include "Dynamics.h"

namespace Engine
{
  class Engine
  {
  public:
    Engine(Simulator* sim) : 
    simulator(sim), t(), eventCount(){}; 
    void equilibrate();
    void simulation();
    
  private:
    Simulator* simulator;
    //pointer to the sampler class
    double t;
    unsigned long long eventCount;

    void handleEvent(Scheduler::Event&, Scheduler::Scheduler&);
    void handleInteraction(Scheduler::Event&, Scheduler::Scheduler&);
    void handleSentinal(Scheduler::Event&, Scheduler::Scheduler&);
    void freeStream(double);
  };
}
