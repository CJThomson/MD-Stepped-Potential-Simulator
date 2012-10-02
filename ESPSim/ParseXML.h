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
 parseXML(SimSet& set, SimProp& prop, SimPotential& pot,std::vector<Particle>& part) :
  settings(set), properties(prop), potential(pot), particles(part) {};
  void parseFile()
  {
    pugi::xml_document config;
    pugi::xml_parse_result result = config.load_file("config.xml");
    if(result)
      {
	parseSettings(config.child("ESPSimConfig").child("SimSettings"));
	parseProperties(config.child("ESPSimConfig").child("SimProperties"));
	parseSampler(config.child("ESPSimConfig").child("SamplerSettings"));
	parseStepper(config.child("ESPSimConfig").child("StepperSettings"));
      }
    //else
    //throw error
  }
 private:
  SimSet& settings;
  SimProp& properties;
  SimPotential& potential;
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
  void parseStepper(pugi::xml_node step)
  {
    pugi::xml_node childNode = step.child("Continuous");
    potential.setContPotential(childNode.attribute("type").as_string(),
			       childNode.attribute("r_cutoff").as_double(),
			       childNode.attribute("epsilon").as_double(),
			       childNode.attribute("sigma").as_double());
    childNode = step.child("Discrete");
    const char* pos = childNode.attribute("position").as_string();
    if(pos == "EvenEnergy")
      potential.setStepPositions(pos, childNode.attribute("energyInterval").as_double());
    else
      potential.setStepPositions(pos, childNode.attribute("noSteps").as_uint());
    potential.setStepEnergies(childNode.attribute("energy").as_string());
    const char* core = childNode.attribute("core").as_string();
    if(core == "Sigma")
      potential.setStepCore(core, childNode.attribute("noSigma").as_uint());
    else if(core == "Manual")
      potential.setStepCore(core, childNode.attribute("corePosition").as_double());
    else
      potential.setStepCore(core);
    potential.setStepPotential(childNode.attribute("potential").as_string());
  }
  void parseProperties(pugi::xml_node prop)
  {
    properties.setT(prop.child("Temperature").attribute("value").as_double());
    properties.setDensity(prop.child("Density").attribute("value").as_double());
    properties.setN(prop.child("NumberOfParticles").attribute("value").as_uint());
  }
};
