#pragma once

#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

#include "Simulator.h"
#include "Particle.h"
#include "Stepmap.h"
#include "Vector3.h"
#include "NeighbourList.h"

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
    inline const Event operator=(const Event &otherevent) 
      {
	collisionTime = otherevent.getCollisionTime();
	particle1 = otherevent.getP1();
	particle2 = otherevent.getP2();
	p2coll = otherevent.getP2Coll();
	eventType = otherevent.getEventType(); 
	return *this;
      }
    inline bool operator<(const Event& otherevent) const  
    { return collisionTime < otherevent.collisionTime; }

    inline bool operator>(const Event& otherevent) const
    { return collisionTime > otherevent.collisionTime; }
	     
    inline double getCollisionTime () const { return collisionTime; }
    inline unsigned int getP1() const { return particle1; }
    inline unsigned int getP2() const { return particle2; }
    inline unsigned long long getP2Coll() const { return p2coll; }
    inline EventType getEventType () const { return eventType; }
  
    inline void setCollisionTime(double value) {collisionTime = value;}
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
  Scheduler(Simulator* sim, boost::shared_ptr<NL::NL> nlist) : 
    simulator(sim), nl(nlist){};
    Event getNextEvent();
    void initialise ();
    void regenerate(double, unsigned long long);
    void update (double, unsigned int); //update for one particles
    void update (double, unsigned int, unsigned int); //update for two particles
    void getThermoEvent(double, unsigned long long);
  private:
    size_t thermoPoint;
    Simulator* simulator;
    std::vector<Event> masterEL;
    Event getMinTime (double, unsigned int);
    double getSentinal (unsigned int);
    double getInteractionTime (unsigned int, unsigned int, Event::EventType&);
    boost::shared_ptr<NL::NL> nl;
    //Schduler::Event getRDF (RDF&);
  };



}
