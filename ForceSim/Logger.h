#pragma once
#include <fstream> //allows output to files
#include <vector> //allows use of vector structures
#include <sstream> //allows use of ostreamstreams
#include "Declares.h" //includes declares for main program
#include "Particle.h" //allows use of particle class
class Logger
{
 public:
  //list of outputs
  typedef enum {
    LOCATIONS,
    BOLTZMANN_H,
    POVRAY
  } OutputType;

  int initialise (OutputType type)
  {
    switch (type)
      {
      case LOCATIONS:
	locLog.open("locLog.dat");
	return 0;
      case BOLTZMANN_H:
	boltzHLog.open("boltzHLog.dat");
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
	if(locLog.is_open())
	  locLog.close();
	return 0;
      case BOLTZMANN_H:
	if(boltzHLog.is_open())
	  boltzHLog.close();
	return 0;
      case POVRAY:
      default:
	return 1;
      };
  }
  void write_BoltzH(double time, double Hvalue)
  {
    if(boltzHLog.is_open())
      boltzHLog << time << "\t" << Hvalue << std::endl;
  }
  void write_Diff (std::vector<Diffusion>& diffusion) //write output for coefficient of diffusion with time
  {
    std::ofstream diffLog;
    diffLog.open ("diffLog.dat"); //open file
    for(int i = 0; i < diffusion.size(); ++i)
      diffLog << diffusion[i].time << "\t" << diffusion[i].coDiff << std::endl; //output all diffusion coeffcients
    diffLog.close(); //close the file
  }
  
  void write_RadDist(double rdfd[], int index, double deltaR,
		      double rho, double T)
  {
    std::ofstream grLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "grLog " << rho << " - " << T << "-forcesim.dat"; 
    outfileName = fileName.str();
    grLog.open(outfileName.c_str());
    for(int i = 0; i < index; ++i)
      grLog << i * deltaR << "\t" << rdfd[i] << "\t" <<  std::endl;
    grLog.close();
  }

  void write_Results(std::vector<Results>& results, double rho, double T)
  {
    std::ofstream resultLog;
    std::ostringstream fileName;
    std::string outfileName;
    fileName << "Results " << rho << " - " << T << "-forcesim.dat"; 
    outfileName = fileName.str();
    resultLog.open(outfileName.c_str());
    for(std::vector<Results>::iterator result = results.begin(); result != results.end(); ++result)
      {
	
	resultLog << result->temperature << "\t"
		  << result->pressure << "\t"
		  << result->potential<< std::endl;
      }
    resultLog.close();
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

  void write_Location (std::vector<CParticle>& particles, double sysTime, CVector3 systemSize)
  {
    if(locLog.is_open())
      {
	locLog << sysTime; //open line with system time
	for(int i = 0; i < particles.size(); ++i)
	  for(int j = 0; j < 3; ++j)
	    locLog << "\t" << particles[i].r[j] - lrint(particles[i].r[j] / systemSize[j])*systemSize[j]; //write each component of each particle's position
	locLog << "\t" << std::endl; //take a new line
      }
  }

 private:
  std::ofstream locLog;
  std::ofstream boltzHLog;
};
