#pragma once
#include <vector>
#include <algorithm>
#include <cmath>

#include "Simulator.h"
#include "Particle.h"
#include "Vector3.h"
namespace NL
{
  class NL
  {
  public:
    virtual void initialise(Simulator* sim) = 0;
    virtual void genNL(unsigned int cell) = 0;
    virtual double getCellTime(unsigned int p) = 0;
    virtual std::vector<unsigned int>& getNeighbours() = 0;
    virtual void moveParticle(unsigned int p) = 0;
  protected:
    NL() {};
  };

  class NL_None: public NL
  {
  public:
    NL_None() {}
    virtual void initialise(Simulator* sim)
    {
      simulator = sim;
      for(std::vector<Particle>::const_iterator p = simulator->getParticles().begin();
	  p != simulator->getParticles().end(); ++p)
	{
	  nl.push_back(p->getID()); //get a vector of every particle ID
	}
    }
    virtual void genNL(unsigned int cell) {}
    virtual double getCellTime(unsigned int p) { return HUGE_VAL;}
    virtual void moveParticle(unsigned int p) {}
    virtual std::vector<unsigned int>& getNeighbours()
      {
	return nl;
      }
  private:
    Simulator* simulator;
    std::vector<unsigned int> nl;

  };

  class NL_Simple: public NL
  {
  public:
  NL_Simple(double rcut) : rmax(rcut) { }
    virtual void initialise(Simulator* sim)
    {
      simulator = sim;
      noCells = floor(simulator->getSysLength() / rmax);
      noCells2 = noCells * noCells;
      cellSize = simulator->getSysLength() / noCells;
      nl.resize(noCells * noCells2);
      generateNeighbourCells();
      for(std::vector<Particle>::const_iterator p = simulator->getParticles().begin();
	  p != simulator->getParticles().end(); ++p)
	{
	  sortParticle(p->getID());
	}
      
    }
    virtual double getCellTime(unsigned int p);
    virtual void moveParticle(unsigned int p) 
    {
      unsigned int cellID = simulator->getParticles()[p].getCell();
      std::vector<unsigned int>::iterator pos = std::find(nl[cellID].begin(),nl[cellID].end(), p);
      if(pos != nl[cellID].end())
	nl[cellID].erase(pos);
      nl[simulator->getParticles()[p].getNextCell()].push_back(p);
      simulator->setParticles()[p].moveCell();
    }
    virtual void genNL(unsigned int cell)
    {
      tempNL.clear();

      for(size_t i(0); i < nc[cell].size(); ++i)
	{
	  /*std::cerr << "i = " << i << " - ";
	  std::cerr << tempNL.size() << " - "; 
	  std::cerr << nl.size() << " - ";
	  std::cerr << nl[i].size() << std::endl;*/
	  tempNL.reserve(tempNL.size() + nl[i].size());
	  tempNL.insert(tempNL.end(), nl[i].begin(), nl[i].end());
	}
      
    }
    virtual std::vector<unsigned int>& getNeighbours()
      {
 	return tempNL;
      }
  private:
    Simulator* simulator;
    std::vector<unsigned int> tempNL;
    std::vector<std::vector<unsigned int> > nl;
    std::vector<std::vector<unsigned int> > nc;
    double rmax;
    double cellSize;
    unsigned int noCells;
    unsigned int noCells2;

    PBCVector<int> conv2Vec(unsigned int) const ;
    unsigned int conv2ID(const PBCVector<int>&) const ;
    void generateNeighbourCells();
    void sortParticle(unsigned int);
  };
}
