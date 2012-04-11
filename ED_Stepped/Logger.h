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

  void write_RadDist (double rdfd[],double rdfc[], int index, double deltaR,
		      double rho, double T)
  {
    std::ofstream grLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "grLog " << rho << " - " << T << "-stepped.dat"; 
    outfileName = fileName.str();
    grLog.open(outfileName.c_str());
    for(int i = 0; i < index; ++i)
      grLog << i * deltaR << "\t" << rdfd[i] << "\t" << rdfc[i] << std::endl;
    grLog.close();
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

  void write_Results (std::vector<Results>& results, double rho, double T)
  {
    std::ofstream resultLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "Results " << rho << " - " << T << "-stepped.dat"; 
    outfileName = fileName.str();
    resultLog.open(outfileName.c_str());
    Results avgResults;
    for(vector<Results>::iterator result = results.begin(); result != results.end(); ++result)
      {
	
	resultLog << result->temperature << "\t"
		  << result->pressure_d << "\t"
		  << result->pressure_c << "\t"
		  << result->potential_d << "\t"
		  << result->potential_c<< endl;
	avgResults += *result;
      }
    avgResults *= 1.0 / results.size();
    resultLog << avgResults.temperature << "\t"
	      << avgResults.pressure_d << "\t"
	      << avgResults.pressure_c << "\t"
	      << avgResults.potential_d << "\t"
	      << avgResults.potential_c << endl;
    resultLog.close();
  }
  void write_Location (std::vector<CParticle>& particles, double sysTime, CVector3 systemSize)
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
