#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "pugixml/pugixml.hpp"
#include "Thermostat/Thermostat.h"
#include "Settings.h"
#include "Sampler.h"

namespace Logger
{
  class Logger
  {
  public:
    void init_Results(bool collisionStats, bool RDFStats)
    {
      runCount = 0;
      results.append_child("ESPSimResults");
      if(collisionStats)
	collStat.append_child("ESPSimCollisionStats");
      if(RDFStats)
	RDFStat.append_child("ESPSimRDFStats");
    }
    void update_Results(Sampler::Sampler& sampler, SimProp properties, bool collisionStats, bool RDFStats)
    {
      ++runCount;
      pugi::xml_node node_Main = results.child("ESPSimResults");
      std::stringstream ss;
      ss << "run" << runCount;
      pugi::xml_node node_Run = node_Main.append_child(ss.str().c_str() );
 
      pugi::xml_node node_Time = node_Run.append_child("Time");
      node_Time.append_attribute("StartTime") = sampler.getStartTime();
      node_Time.append_attribute("EndTime") = sampler.getEndTime();
      node_Time.append_attribute("Events/sec") = sampler.getEventPS();
      node_Time.append_attribute("SimTime/sec") = sampler.getTimePS();
     
      pugi::xml_node node_Length = node_Run.append_child("Length");
      node_Length.append_attribute("Time") = sampler.getTime();
      node_Length.append_attribute("EventCount") = sampler.getEventCount();
      node_Length.append_attribute("MeanFreeTime") = sampler.getMeanFreeTime();

      pugi::xml_node node_Temp = node_Run.append_child("Temperature");
      node_Temp.append_attribute("Mean") = sampler.getMeanT();
      node_Temp.append_attribute("MeanSqr") = sampler.getMeanSqrT();
      node_Temp.append_attribute("StdDev") = sqrt(sampler.getMeanSqrT() 
						  - pow(sampler.getMeanT(), 2));

      pugi::xml_node node_U = node_Run.append_child("Potential");
      node_U.append_attribute("Mean") = sampler.getMeanU();
      node_U.append_attribute("MeanSqr") = sampler.getMeanSqrU();
      node_U.append_attribute("StdDev") = sqrt(sampler.getMeanSqrU() 
					       - pow(sampler.getMeanU(), 2));
      pugi::xml_node node_P = node_Run.append_child("Pressure");
      node_P.append_attribute("Mean") = sampler.getP();
 
      if(collisionStats)
	write_Collisions(sampler);
      if(RDFStats)
	write_RDF(sampler, properties);
    }
    void write_Results(bool collisionStats, bool RDFStats)
    {
      saveXML(results,"Results.xml");
      if(collisionStats)
	saveXML(collStat,"CollisionStatistics.xml");
      if(RDFStats)
	saveXML(RDFStat, "RDFStatistics.xml");
    }
    void write_Collisions(Sampler::Sampler& sampler)
    {
      std::stringstream ss;
      ss << "run" << runCount;
      pugi::xml_node node_Main = collStat.child("ESPSimCollisionStats");
      pugi::xml_node node_Run = node_Main.append_child(ss.str().c_str() );
      for(std::vector<Sampler::CollisionCount>::const_iterator i = sampler.getCollCount().begin(); 
	  i != sampler.getCollCount().end(); ++i)
	{
	  pugi::xml_node node_step = node_Run.append_child("Step");
	  node_step.append_attribute("ID") = i->getID();
	  pugi::xml_node node_stepProp = node_step.append_child("StepProperties");
	  node_stepProp.append_attribute("R") = i->getR();
	  node_stepProp.append_attribute("U") = i->getU();
 	  pugi::xml_node node_collCapt = node_step.append_child("Captures");
	  node_collCapt.append_attribute("In") = i->getCount(2,true);
	  node_collCapt.append_attribute("Out") = i->getCount(2,false);
	  pugi::xml_node node_collBounce = node_step.append_child("Bounces");
	  node_collBounce.append_attribute("In") = i->getCount(4,true);
	  node_collBounce.append_attribute("Out") = i->getCount(4,false);
	}
    }
    void write_RDF(Sampler::Sampler& sampler, SimProp properties)
    {
      std::stringstream ss;
      ss << "run" << runCount;
      pugi::xml_node node_Main = RDFStat.child("ESPSimRDFStats");
      pugi::xml_node node_Run = node_Main.append_child(ss.str().c_str() );
      std::vector<double> RDF_data = sampler.calcRDF(properties.getN(), properties.getDensity());
      double binWidth = sampler.getBinWidth();
      for(size_t i = 0; i < RDF_data.size(); ++i)
	{
	  pugi::xml_node node_Bin = node_Run.append_child("Bin");
	  node_Bin.append_attribute("R") = i * binWidth;
	  node_Bin.append_attribute("RDF") = RDF_data[i];
	}

    }
    void write_outConfig(SimSet& simSettings, SimProp& simProperties, SimPotential& simPot)
      {
	pugi::xml_document outConfig;
	pugi::xml_node node_Main = outConfig.append_child("ESPSimConfig");
	pugi::xml_node node_SimProp = node_Main.append_child("SimProperties");
	
	pugi::xml_node node_Temp = node_SimProp.append_child("Temperature");
	node_Temp.append_attribute("value") = simProperties.getT();

	pugi::xml_node node_Density = node_SimProp.append_child("Density");
	node_Density.append_attribute("value") = simProperties.getDensity();

	pugi::xml_node node_NPart = node_SimProp.append_child("NumberOfParticles");
	node_NPart.append_attribute("value") = simProperties.getN();
	
	pugi::xml_node node_SimSet = node_Main.append_child("SimSettings");

	pugi::xml_node node_EqLength = node_SimSet.append_child("EqLength");
	if(simSettings.isTime(true))
	  node_EqLength.append_attribute("time") = simSettings.getEQTime();
	else
	  node_EqLength.append_attribute("event") = (unsigned int)simSettings.getEQEvent();

	pugi::xml_node node_RunLength = node_SimSet.append_child("RunLength");
	if(simSettings.isTime(false))
	  node_RunLength.append_attribute("time") = simSettings.getRunTime();
	else
	  node_RunLength.append_attribute("event") = (unsigned int)simSettings.getRunEvent();
	pugi::xml_node node_Thermostat = node_SimSet.append_child("Thermostat");
	node_Thermostat.append_attribute("type") = simSettings.getThermostat()->getType();
	node_Thermostat.append_attribute("AutoUpdate") = simSettings.getThermoControl();
	node_Thermostat.append_attribute("ThermoFreq") = simSettings.getThermoFreq();
	pugi::xml_node node_NL = node_SimSet.append_child("NL");
	node_NL.append_attribute("type") = simSettings.getNLType();

	pugi::xml_node node_Sampler = node_Main.append_child("SamplerSettings");
	
	pugi::xml_node node_CollCount = node_Sampler.append_child("CollisionCounts");
	node_CollCount.append_attribute("active") = simSettings.getSampleColl();

	pugi::xml_node node_Potential = node_Main.append_child("StepperSettings");
	
	pugi::xml_node node_Cont = node_Potential.append_child("Continuous");
	node_Cont.append_attribute("type") = simPot.getContPotential();
	node_Cont.append_attribute("epsilon") = simPot.getEpsilon();
	node_Cont.append_attribute("sigma") = simPot.getSigma();
	node_Cont.append_attribute("r_cutoff") = simPot.getRCutOff();

	pugi::xml_node node_DisCont = node_Potential.append_child("Discrete");
	node_DisCont.append_attribute("position") = simPot.getStepPositions();
	node_DisCont.append_attribute("energy") = simPot.getStepEnergies();
	node_DisCont.append_attribute("core") = simPot.getStepCore();
	if(strcmp(simPot.getStepPositions(), "EvenEnergy") == 0)
	  node_DisCont.append_attribute("energyInterval") = simPot.getEnergyInterval();
	else
	  node_DisCont.append_attribute("noSteps") = simPot.getNoStep();
	
	saveXML(outConfig, "Out.Config.xml");
      }
  private:
    pugi::xml_document results; //results document
    pugi::xml_document collStat; //collision statistics document
    pugi::xml_document RDFStat; //RDF statistics document
    unsigned int runCount;
    void saveXML(pugi::xml_document& doc, const char* name)
    {
      doc.save_file(name);
    }
  };

}
