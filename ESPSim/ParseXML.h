#pragma once

#include <vector>
#include <boost/shared_ptr.hpp>
#include "pugixml/pugixml.hpp"
#include "Particle.h"
#include "Settings.h"
#include "Thermostat/Thermostat.h"
class parseXML
{
 public: 
 parseXML(SimSet& set, SimProp& prop, std::vector<Particle>& part) :
  settings(set), properties(prop), particles(part) {};
  void parseFile()
  {
    pugi::xml_document config;
    pugi::xml_parse_result result = config.load_file("config.xml");
    if(result)
      {
	parseSettings(config.child("ESPSimConfig").child("SimSettings"));
	parseProperties(config.child("ESPSimConfig").child("SimProperties"));
	parseSampler(config.child("ESPSimConfig").child("SamplerSettings"));
      }
    //else
    //throw error
  }
 private:
  SimSet& settings;
  SimProp& properties;
  std::vector<Particle> particles;

  void parseSettings(pugi::xml_node set)
  {
    //Equilibrium/Production Lengths
    settings.setRunEvent(set.child("RunLength").attribute("event").as_uint());
    settings.setRunTime(set.child("RunLength").attribute("time").as_double());
    settings.setEQEvent(set.child("EqLength").attribute("event").as_uint());
    settings.setEQTime(set.child("EqLength").attribute("time").as_double());
    //Thermostat
    settings.setThermoControl(set.child("Thermostat").attribute("autoUpdate").as_bool());
    settings.setThermoFreq(set.child("Thermostat").attribute("thermoFreq").as_double());
    settings.setThermoType(set.child("Thermostat").attribute("type").as_string());

  }
  void parseSampler(pugi::xml_node samp)
  {
    settings.setSampleColl(samp.child("CollisionCounts").attribute("active").as_bool());
  }
  void parseProperties(pugi::xml_node prop)
  {
    properties.setT(prop.child("Temperature").attribute("value").as_double());
    properties.setDensity(prop.child("Density").attribute("value").as_double());
    properties.setN(prop.child("NumberOfParticles").attribute("value").as_uint());
  }
};
