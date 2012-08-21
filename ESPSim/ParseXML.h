#pragma once

#include <vector>

#include "pugixml/pugixml.hpp"
#include "Particle.h"
#include "Settings.h"
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
    settings.setRunEvent(set.child("RunLength").attribute("event").as_uint());
    settings.setRunTime(set.child("RunLength").attribute("time").as_double());
    settings.setEQEvent(set.child("EqLength").attribute("event").as_uint());
    settings.setEQTime(set.child("EqLength").attribute("time").as_double());
  }

  void parseProperties(pugi::xml_node prop)
  {
    properties.setT(prop.child("Temperature").attribute("value").as_double());
    properties.setDensity(prop.child("Density").attribute("value").as_double());
    properties.setN(prop.child("NumberOfParticles").attribute("value").as_uint());
  }
};
