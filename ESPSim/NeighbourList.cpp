#include "NeighbourList.h"

namespace NL
{
  
  void NL_Simple::generateNeighbourCells()
    {
      nc.clear(); 
      nc.resize(noCells * noCells2);

      for(size_t i (0); i < nc.size(); ++i) 
	{

	  PBCVector<int> cell = conv2Vec(i);
	  
	  for(int dx = -1; dx < 2; ++dx) 
	    for(int dy = -1; dy < 2; ++dy)
	      for(int dz = -1; dz < 2; ++dz)
		{
		  PBCVector<int> cell2(noCells, false, 
						Vector3<int> (dx, dy, dz));
		  unsigned int cellID = conv2ID(cell + cell2);
		  if(std::find(nc[i].begin(), nc[i].end(), cellID) == nc[i].end())
		    nc[i].push_back(cellID);
		}
	}
    }
  
  void NL_Simple::sortParticle(unsigned int p)
   {
     PBCVector<int> cell(noCells, false, Vector3<int>());
     PBCVector<double> r(simulator->getSysLength(), true, 
			 simulator->getParticles()[p].getR());
     for(size_t i(0); i < 3; ++i) 
       cell[i] = floor((r[i] + 0.5 * simulator->getSysLength()) / cellSize); 
     unsigned int cellID = conv2ID(cell);
     nl[cellID].push_back(p);
     simulator->setParticles()[p].setCell(cellID);
   }
  double NL_Simple::getCellTime(unsigned int p)
  {
    PBCVector<double> location(simulator->getSysLength(), true, 
			       simulator->getParticles()[p].getR()); 
    Vector3<double> t_min(HUGE_VAL, HUGE_VAL, HUGE_VAL); 
    Vector3<double> boundary;
    PBCVector<int> cell(conv2Vec(simulator->getParticles()[p].getCell())); 
    PBCVector<int> direction(noCells, false, Vector3<int>());
    for(size_t i(0); i < 3; ++i) 
      {
	if(simulator->getParticles()[p].getV()[i] != 0)
	  {
	    if(simulator->getParticles()[p].getV()[i] < 0) //if particle is going backwards
	      {
		boundary[i] = cell[i] * cellSize - 0.5 * simulator->getSysLength(); 
		direction[i] = -1;
	      }
	    else
	      {
		boundary[i] = (cell[i] + 1) * cellSize - 0.5 * simulator->getSysLength();
		direction[i] = 1;
	      }
	  }

	double distance = boundary[i] - location[i]; 
	t_min[i] = fabs(distance / simulator->getParticles()[p].getV()[i]);     
	/*std::cerr << "p: " << p
		  << " b - l = d: " << boundary[i] 
		  << " - " << location[i]
		  << " = " << distance 
		  << " v: " << simulator->getParticles()[p].getV()[i]
		  << " t: " << t_min[i] << std::endl;*/
      }
    
    double tMin =  std::min(t_min[0], std::min(t_min[1], t_min[2])); 
    for(size_t i(0); i < 3; ++i)
      if(tMin != t_min[i])
	direction[i] = 0;

    PBCVector<int> sum (noCells, false, cell + direction);
    /*std::cerr << "sum = (" << sum[0] << ", " << sum[1] 
      << ", " << sum[2] << ")" << std::endl;*/
    unsigned int nextCell = conv2ID(sum);
    if(nextCell > noCells * noCells2)
      {
	std::cerr << "Error: Invalid next cell" << std::endl;
	exit(34);
      }
    simulator->setParticles()[p].setNextCell(nextCell);
    return tMin; //return the minimum time
  }

  PBCVector<int> NL_Simple::conv2Vec(unsigned int ID) const
  {
    int x = int(ID / noCells2);
    ID -= x * noCells2;
    int y = int(ID / noCells);
    ID -= y * noCells;
    int z = ID;
    return PBCVector<int>(noCells, false, Vector3<int>(x, y,z));
  }
  
  unsigned int NL_Simple::conv2ID(const PBCVector< int>& vec) const
  {
    return vec[0] * noCells2 + vec[1] * noCells + vec[2];
  }
}
