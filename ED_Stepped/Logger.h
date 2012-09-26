#pragma once
#include <fstream> //allows output to files
#include <vector> //allows use of vector structures
#include <string>
#include "Declares.h" //includes declares for main program
#include "Particle.h" //allows use of particle class
class Logger
{
 public:
  //list of outputs
  typedef enum {
    LOCATIONS,
    OUTPUTLOG,
    POVRAY
  } OutputType;

  int initialise (OutputType type)
  {
    switch (type)
      {
      case LOCATIONS:
	locLog.open("locLog.dat");
    startLoc = false;
	return 0;
	case OUTPUTLOG:
    outLog.open("outLog.dat");
  return 0;
      case POVRAY:
      default:
	return 1;
      };
  }
  int terminate (OutputType type)
  {
    switch (type)
      {
      case LOCATIONS:
	locLog.close();
	return 0;
    case OUTPUTLOG:
      if(outLog.is_open()) { outLog.close();}
      case POVRAY:
      default:
	return 1;
      };
  }

  void write_Diff (std::vector<Diffusion>& diffusion) //write output for coefficient of diffusion with time
  {
    std::ofstream diffLog;
    diffLog.open ("diffLog.dat"); //open file
    for(int i = 0; i < diffusion.size(); ++i)
      diffLog << diffusion[i].time << "\t" << diffusion[i].coDiff << std::endl; //output all diffusion coeffcients
    diffLog.close(); //close the file
  }

  void write_RadDist (double rdfd[], int index, double deltaR,
		      double rho, double T, double n, double stepInt = -1)
  {
    std::ofstream grLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "grLog " << rho << " - " << T << " - " << n << " - " << stepInt <<"-stepped.dat";
    outfileName = fileName.str();
    grLog.open(outfileName.c_str());
    for(int i = 0; i < index; ++i)
      grLog << i * deltaR << "\t" << rdfd[i] << std::endl;
    grLog.close();
  }
  void writeScriptInput(std::vector<Steps>& steps, std::vector<CollisionCount>& eCount, unsigned int stepNo)
  {
    std::ofstream scriptInput;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << stepNo << "Steps.dat";
    outfileName = fileName.str();
    scriptInput.open(outfileName.c_str());
    for(size_t i(0); i < steps.size(); ++i)
      {
	scriptInput << steps[i].step_radius << "\t" << steps[i].step_energy  << "\t"
		    << eCount[i].in_bounce + eCount[i].out_bounce << std::endl;
      }
    scriptInput.close();
  }
  void write_Steps(std::vector<Steps>& steps, double T, double rho, double n, 
		   Stepper::StepHeight type, double stepInt= -1)
  {
    std::ofstream stepLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "stepLog " << rho << " - " << T << " - " << n << " - " << type 
	     << " - " << stepInt << "-stepped.dat";
    outfileName = fileName.str();
    stepLog.open(outfileName.c_str());
    for(std::vector<Steps>::iterator i = steps.begin(); i != steps.end(); ++i)
      {
	stepLog << i->step_radius << "\t" << i->step_energy << std::endl;
      }
    stepLog.close();
  }
  void write_contRDF(std::vector<std::pair<double, double> >& contRDF, double rho, 
		     double T, double n, double stepInt = -1)
  {
    std::ofstream grLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "contRDFLog " << rho << " - " << T << " - " << n << " - " << stepInt <<"-stepped.dat";
    outfileName = fileName.str();
    grLog.open(outfileName.c_str());
    for(std::vector<std::pair<double, double> >::iterator i = contRDF.begin(); i != contRDF.end(); ++i)
      {
	grLog << i->first << "\t" << i->second << std::endl;
      }
    grLog.close();
  }
  void write_ICF(std::vector<std::pair<double, double> >& icf, double rho, 
		 double T, double n, double stepInt = -1)
  {
    std::ofstream icfLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "ICFLog " << rho << " - " << T << " - " << n << " - " << stepInt <<"-stepped.dat";
    outfileName = fileName.str();
    icfLog.open(outfileName.c_str());
    for(std::vector<std::pair<double, double> >::iterator i = icf.begin(); i != icf.end(); ++i)
      {
	icfLog << i->first << "\t" << i->second << std::endl;
      }
    icfLog.close();
  }

  void write_Init (std::vector<CParticle>& particles)
  {
    std::ofstream initLog;
    initLog.open("init.dat");
    for(int i = 0; i < particles.size(); ++i)
      {
	for(int j = 0; j < 3; ++j)
	  initLog << particles[i].r[j] << "\t";
	for(int j = 0; j < 3; ++j)
	  initLog << particles[i].v[j] << "\t";
	initLog << std::endl;
      }
    initLog.close();
  }

  void write_CollCount(std::vector<CollisionCount>& eCount, double rho, 
		       double T, double n, double stepInt = -1)
  {
    std::ofstream collLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "collCountLog " << rho << " - " << T << " - " << n << " - " << stepInt <<"-stepped.dat";
    outfileName = fileName.str();
    collLog.open(outfileName.c_str());
    for(std::vector<CollisionCount >::iterator i = eCount.begin(); 
	i != eCount.end(); ++i)
      {
	collLog << i->in_capture << "\t" << i->in_bounce << "\t"
		<< i->out_release << "\t" << i->out_bounce << std::endl;
      }
    collLog.close();
  }

  void write_Results(std::vector<Results>& results, double rho, double T, int n, int noEvents, double stepInt = -1)
  {
    std::ofstream resultLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "Results " << rho << " - " << T << "-" << n << " - " << stepInt << "-stepped.dat";
    outfileName = fileName.str();
    resultLog.open(outfileName.c_str());
    for(std::vector<Results>::iterator result = results.begin(); result != results.end(); ++result)
      {
	resultLog << result->temperature << "\t"
		  << result->noEvents << "\t"
		  << result->pressure_d << "\t"
		  << result->pressure_c << "\t"
		  << result->potential_d << "\t"
		  << result->potential_c << std::endl;
      }

    resultLog.close();
  }
  void write_Location (std::vector<CParticle>& particles, double sysTime, CVector3<double> systemSize)
  {
    if(locLog.is_open())
      {
        if(!startLoc)
        {
                  locLog << systemSize[0] << "\t" << systemSize[1] << "\t" <<systemSize[2] << "\t" << std::endl;
                  startLoc = true;
        }

	locLog << sysTime; //open line with system time
	for(int i = 0; i < particles.size(); ++i)
	  for(int j = 0; j < 3; ++j)
	    locLog << "\t" << particles[i].r[j]; //- lrint(particles[i].r[j] / systemSize[j])*systemSize[j]; //write each component of each particle's position
	locLog << "\t" << std::endl; //take a new line
      }
  }
  std::ofstream outLog;
 private:
  std::ofstream locLog;
  bool startLoc;
  
};
