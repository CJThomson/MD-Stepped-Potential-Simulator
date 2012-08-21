#pragma once

#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
#include "Simulator.h"
#include "Particle.h"
#include "Stepmap.h"
#include "Vector3.h"

typedef std::vector<Particle>::iterator it_p;
namespace Scheduler
{

  class Event
  {
  public:	
    typedef enum 
    {
      IP_OUT = 0,
      IP_IN = 1,
      SENTINAL = 2,
      NEIGHBOURCELL = 3,
      THERMOSTAT = 4,
      RDF = 5,
      NONE = 6
    }EventType;

  Event(double dt, int p1, int p2, int coll, EventType _type):
    collisionTime(dt),
      particle1(p1),
      particle2(p2),
      p2coll(coll),
      eventType(_type)
      {}

  Event():
    collisionTime(HUGE_VAL),
      particle1(-1),
      particle2(-1),
      p2coll(0),
      eventType(NONE)
      {}

    bool operator<(const Event& otherevent) const  
    { return collisionTime < otherevent.collisionTime; }

    bool operator>(const Event& otherevent) const
    { return collisionTime > otherevent.collisionTime; }
	     
    double getCollisionTime () {return collisionTime;}
    unsigned int getP1() {return particle1;}
    unsigned int getP2() {return particle2;}
    unsigned long long getP2Coll() { return p2coll; }
    EventType getEventType () { return eventType; }
  
    void setCollisionTime(double value) {collisionTime = value;}
  private:
    double collisionTime;
    unsigned int particle1;
    unsigned int particle2;
    unsigned long long p2coll;
    EventType eventType;
  };

  class Scheduler
  {
  public: 
    Scheduler(Simulator* sim) : 
    simulator(sim) {};
    Event getNextEvent();
    void initialise ();
    void update (double, unsigned int); //update for one particles
    void update (double, unsigned int, unsigned int); //update for two particles
  private:
    Simulator* simulator;
    std::vector<Event> masterEL;
    Event getMinTime (double, unsigned int);
    double getSentinal (unsigned int);
    double getInteractionTime (unsigned int, unsigned int, Event::EventType&);
    //Schduler::Event getThermostat (Thermostat&);
    //Schduler::Event getRDF (RDF&);
  };



}
