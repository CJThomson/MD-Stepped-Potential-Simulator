#pragma once

#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "Particle.h"
#include "Simulator.h"
#include "Scheduler.h"
#include "Dynamics.h"
#include "Logger.h"
namespace Engine
{
  class Engine
  {
  public:
    Engine(Simulator* sim) : 
    simulator(sim), t(), eventCount(){}; 
    void equilibrate();
    void productionRun();
    void simulation();
    
  private:
    Simulator* simulator;
    //pointer to the sampler class
    double t;
    unsigned long long eventCount;
    bool equilibration;
    void handleEvent(Scheduler::Event&, Scheduler::Scheduler&, Sampler::Sampler& );
    void handleInteraction(Scheduler::Event&, Scheduler::Scheduler&, Sampler::Sampler& );
    void handleSentinal(Scheduler::Event&, Scheduler::Scheduler&, Sampler::Sampler&);
    void freeStream(double, Sampler::Sampler&);
  };
}
