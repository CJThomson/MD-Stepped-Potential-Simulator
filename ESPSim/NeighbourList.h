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
    virtual double getCellTime(unsigned int p) = 0;
    virtual std::vector<unsigned int>& getNeighbourCells(unsigned int cell) = 0;
    virtual std::vector<unsigned int>& getParticles(unsigned int cell) = 0;
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
      nc.push_back(0);
    }
    virtual std::vector<unsigned int>& getNeighbourCells(unsigned int cell)
      {
	return nc;
      }
    virtual std::vector<unsigned int>& getParticles(unsigned int cell)
      {
	return nl;
      }
    virtual double getCellTime(unsigned int p) { return HUGE_VAL;}
    virtual void moveParticle(unsigned int p) {}

  private:
    Simulator* simulator;
    std::vector<unsigned int> nl;
    std::vector<unsigned int> nc;
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
    virtual std::vector<unsigned int>& getNeighbourCells(unsigned int cell)
      {
	return nc[cell];

      }

    virtual std::vector<unsigned int>& getParticles(unsigned int cell)
      {
	return nl[cell];
      }
    virtual void moveParticle(unsigned int p) 
    {
      unsigned int cellID = simulator->getParticles()[p].getCell();
      std::vector<unsigned int>::iterator pos = std::find(nl[cellID].begin(),nl[cellID].end(), p);
      if(pos != nl[cellID].end())
	nl[cellID].erase(pos);
      nl[simulator->getParticles()[p].getNextCell()].push_back(p);
      simulator->setParticles()[p].moveCell();
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
