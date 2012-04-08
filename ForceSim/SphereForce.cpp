//----Program Includes----
#include "Declares.h" //declarations for this program

//Physical Properties:
const double density = 0.776;
const double temperature = 0.85; //temperature of the system

//Simulation:
int numberParticles = 864; //number of particles
const int simTime = 60000; //length of the Simulation
const double dt = 0.005; //length of ticme interval
const double length = pow(numberParticles / density, 1.0 / 3.0);
const double r_cut = 3.0;
const CVector3 systemSize(length, length, length); //size of the system
const bool initFile = false; //use an init file
const bool overwriteInit = false; //create a new initfile
const int NL_update = 10;

const double ptail = -16 * M_PI * density * density / (3.0 * pow(r_cut, 3)) * (1.0 - 2.0 / (3.0 * pow(r_cut,6))); 
const double utail = -8 * M_PI * density / (3.0 * pow(r_cut, 3)) * (1.0 - 1.0 / (3.0 * pow(r_cut,6))); 

//Thermostat
bool thermostat = true; //use a thermostat
const double thermoFreq = 0.05; //update frequency of thermostat
const int thermoOff = 10000;

//Reduced Unit Definitions:
const double mass = 1; //mass of a particle
const double radius = 0.5; //radius of a particle (set for diameter = 1)
const double epsilon = 1; //minimum energy of Lennard Jones Potential
const double sigma = 1; //distance for Lennard Jones root

//Logging:
const int out_interval = 20; //frequency of output to file
const int diff_interval = 10;
const int rdf_interval = 10;
const int sample_interval = 10;
const bool writeLoc = false;

//Measuring Properties:
const int startSampling = 20000; //number of readings to take
double readingTime = 0;
const int noBins = 500; //number of radial bins
const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
bool calcP = true;

double gVal[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over time

//----Time Averages----
double TA_gVal[noBins];
double TA_Temp;
double TA_U;
double TA_Pressure;

using namespace std;

int main()
{
  //variable declarations
  Logger log; //create instance of the logger class
  vector<CParticle> particles;
  int noReadings = 0;
  cout << "Initialising random number generator..." << endl;
  CRandom RNG;

  //Initialise the simulation
  if(initFile)
    initFromFile(particles); //initialise the system from file
  else
    initialise(particles, RNG); // initialise the system
  numberParticles = particles.size();
  cout << "Initialising " << numberParticles << " particles with density "
       << density << " is a box of side length " << length << endl;
  vector<int> neighbourList;
  int listPos[particles.size()];

  cout << "Initialising Log files..." << endl;
  if(writeLoc) { log.initialise(Logger::LOCATIONS); } //initialise location logger

  cout << "Initialising neighbour lists..." << endl;
  calcNeighbourList(particles, neighbourList, listPos);
  calcAllForces(particles, neighbourList, listPos); //calculate all forces on system

  if(writeLoc) { log.write_Location(particles, 0, systemSize); } //write initial values to log

  cout << "Starting Simulation..." << endl;
  double t = 0;
  for (size_t timeStep (0); timeStep < simTime; ++timeStep)
    {
      t += dt;
      //move all particles forward
      for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
	{
	  //move particles in space using Velocity Verlet
	  particle->r += dt * (particle->v + particle->a * dt * 0.5);
	  particle->v += dt * particle->a * 0.5; //calculate part of velcity
	  (particle->a).zero(); //zero acceleration
	}

      calcAllForces(particles, neighbourList, listPos); //calculate forces (and accelerations)

      for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
	particle->v += dt * particle->a * 0.5; //calculate remainder of velocity

      //apply Andersen Thermostat
      if(thermostat)
	{
	  double factor = sqrt(temperature);
	  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
	    if(RNG.var_01() < thermoFreq)
	      {
		for(int j = 0; j < 3; ++j)
		  particle->v[j] = RNG.var_normal() * factor;
	      }
	}

      //update particle
      if(timeStep % out_interval == 0)
	log.write_Location(particles, t, systemSize); //write locations to output file

      if(timeStep % NL_update == 0)
	calcNeighbourList(particles, neighbourList, listPos);

      if(timeStep % 100 == 0) //write a screen output every 100 timesteps
	{
	  double pshort = (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3)) + ptail;
      	  cout << fixed << setprecision(2) << (double) timeStep / simTime * 100 << "% complete" << "\t"
	       << " t: " << setprecision(2) << t << "\t"
	       << " T: " << setprecision(4) << calcTemp(particles) << "\t"
	       << " P: " << setprecision(3) << pshort << "\t"
	       << " P + LR: " << setprecision(3) << pshort + ptail << "\t"
	       << " U: " << setprecision(3) << calcPotential(particles,neighbourList, listPos) << "\t"
	       << endl;
	}

      if(timeStep == thermoOff && thermostat)
      {
	thermostat = false;  //turn off thermostat just prior to taking readings
	zeroMomentum(particles); //zero linear momentum in the system

      }
      if(timeStep == startSampling)
	{
	  readingTime = t;
	  for(int i = 0; i < particles.size(); ++i)
	    particles[i].r0 = particles[i].r;
	}

      if(timeStep < startSampling)
	{
	  ++noReadings;
	  if(timeStep % diff_interval == 0)
	    {
	      double D = calcDiff(particles, t); //calculate value of diffusion coefficient
	      coDiff.push_back(Diffusion(D,t)); //add value to vector
	    }

	  if(timeStep % rdf_interval == 0)
	    {
	      calcRadDist(particles); //calculate radial distribution
	      for(int i = 0; i < noBins; ++i)
		TA_gVal[i] += gVal[i]; //calculate running total
	    }

	  if(timeStep % sample_interval == 0)
	    {
	      TA_Pressure += (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3)); //add current pressure to TA
	      TA_Temp += calcTemp(particles); //add current temp to TA
	      TA_U += calcPotential(particles, neighbourList, listPos); //add current potential energy to TA
	    }
	}

    }

  cout << "Simulation Complete" << endl;
  //After Simulation
  double deltaR = maxR / noBins; //width of each shell
  for(size_t i(0); i < noBins; ++i)
    {
      double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
      TA_gVal[i]/=(0.5 * numberParticles * noReadings / rdf_interval * volShell * density);
    }
  cout << "write output files..." << endl;
  log.write_Location(particles, simTime, systemSize);
  log.write_RadDist(TA_gVal,noBins, deltaR); //write radial gdistribution file
  log.write_Diff(coDiff); //write coefficient of diffusion file
  if(overwriteInit)
    log.write_Init(particles);

  //output Time Averaged system properties
  cout << setprecision(5)
       << "Time Averages:- Temp: " << TA_Temp * sample_interval / noReadings
       << " P: " << TA_Pressure * sample_interval / noReadings
       << " P+LR: " << TA_Pressure * sample_interval / noReadings + ptail 
       << " U: " << TA_U * sample_interval / (numberParticles * noReadings) 
       << " U+LR: " << TA_U * sample_interval / (numberParticles * noReadings) + utail
       << endl;
  return 0;
}
								    
void initFromFile (vector<CParticle> &particle)
{
  ifstream initLog ("init.dat");
  string line;
  int strpos = 1;
  while(initLog.good() )
    {
      int counter = 0; //rezero counter
      CVector3 init_r;
      CVector3 init_v;
      getline(initLog, line); //get next line of log
      vector<CParticle> p;
      for(size_t i (0); i < line.size(); ++i)
        {
	  if(line[i] == '\t')
            {
	      switch(counter)
                {
                case 0: //then position x
		  init_r.x = atof(line.substr(0, i).c_str());
		  strpos = i+1;
		  break;
                case 1: //then position y
		  init_r.y = atof(line.substr(strpos,i - strpos).c_str());
		  strpos = i+1;
		  break;
                case 2: //then position z
		  init_r.z = atof(line.substr(strpos,i - strpos).c_str());
		  strpos = i+1;
		  break;
                case 3: //then velocity x
		  init_v.x = atof(line.substr(strpos,i - strpos).c_str());
		  strpos = i+1;
		  break;
                case 4: //then velocity y
		  init_v.y = atof(line.substr(strpos,i - strpos).c_str());
		  strpos = i+1;
		  break;
                case 5: //then z
		  init_v.z = atof(line.substr(strpos,i - strpos).c_str());
		  particle.push_back(CParticle(init_r, init_v, radius, mass, i));
		  break;
                default:
		  exit(1);
		  cout << "error in reading file" << endl;
                }
	      ++counter;
            }
        }
    }
}

void initialise(vector<CParticle>& particle, CRandom& rng)
{
  int n = ceil(pow(numberParticles / 4, (double) 1 / 3)); //find cubic root of number of particles
  double a = length / n;
  int j = 0, x = 0, y = 0, z = 0;
  for(int i = 0; i < numberParticles; ++i)
    {
      CVector3 location;
      switch (j) //set particles in a FCC structure
	{
	case 0:
	  location = CVector3(x * a - length * 0.5,
			      y * a - length * 0.5,
			      z * a - length * 0.5);
	  break;
	case 1:
	  location = CVector3(x * a - length * 0.5,
			      (y + 0.5) * a - length * 0.5,
			      (z + 0.5) * a - length * 0.5);
	  break;
	case 2:
	  location = CVector3((x + 0.5) * a - length * 0.5,
			      y * a - length * 0.5,
			      (z + 0.5) * a - length * 0.5);
	  break;
	case 3:
	  location = CVector3((x + 0.5) * a - length * 0.5,
			      (y + 0.5) * a - length * 0.5,
			      z * a - length * 0.5);
	  j = -1;
	  break;
	}
      ++j;
      if(j==0)
	++x;
  
      if(x >= n)
	{
	  x=0;
	  ++y;
	}
      if(y >= n)
	{
	  y=0;
	  ++z;
	}
      
      CVector3 velocity(rng.var_normal(),
			rng.var_normal(),
			rng.var_normal()); //assign a random unit vector to the velocity
      particle.push_back(CParticle(location, velocity, radius, mass, i));
    }

  correctVelocity(particle);
}

double calcDiff(vector<CParticle>& particle, double time)
{
  double sumDiff = 0;
  for(int i = 0; i < particle.size(); ++i)
    {
      CVector3 distTravelled = particle[i].r - particle[i].r0;
      sumDiff += distTravelled.dotProd(distTravelled);
    }


  return sumDiff / particle.size();
  //return sumDiff / (6 * particle.size() * (time -readingTime));
}
double calcTemp(vector<CParticle> &particle)
{
  double sum = 0;
  for(int i = 0; i < particle.size(); ++i)
    sum += particle[i].v.dotProd(particle[i].v);
  return mass/(3*particle.size())*sum;
}
void calcRadDist(vector<CParticle> &particle)
{
  for(int i = 0; i < noBins; ++i)
    gVal[i] = 0; //rezero array

  for(int i = 0; i < particle.size(); ++i)
    {
      for(int j = i + 1; j < particle.size(); ++j)
	{
	  CVector3 distance = particle[i].r - particle[j].r;
	  applyBC(distance);
	  if(distance.length() < maxR)
	    {
	      int index = floor(distance.length() * noBins/ maxR);
	      ++gVal[index];
	    }
	}
    }
}

CVector3 calcForce(CVector3 rij)
{
  return 24.0 * epsilon / sigma * (2.0 * pow(sigma / rij.length(), 13) - pow(sigma / rij.length(), 7)) * rij.normalise();
}
double calcKinetic(vector<CParticle> &particle)
{
  double kinetic=0;
  for(int i = 0; i < numberParticles; i++)
    kinetic += particle[i].kineticEnergy();
  return kinetic;
}
double calcPotential(vector<CParticle> &particle, vector<int> &NL, int NLpos[])
{
  double potential = 0;
  //loop through each pair of particles
  for(int i = 0; i < particle.size(); ++i)
    {
      int max_length = NL.size();
      if(i+1 < particle.size())
        max_length = NLpos[i+1];
      for(int j = NLpos[i]; j < max_length; ++j)
	{
	  CVector3 distance = particle[i].r - particle[NL[j]].r; //vector between particles
	  applyBC(distance);
	  if(distance.length() <= 6 * radius)
      potential += 4.0 * epsilon*(pow(sigma/distance.length(),12)-pow(sigma/distance.length(),6));
	}
    }
  return potential;
}
void correctVelocity(vector<CParticle> &particle)
{
  zeroMomentum(particle); //zero linear momentum of the system
  double temp = calcTemp(particle); //calculate temperature of system
  double factor = sqrt(temperature/temp); //calculate correction factor to set system to correct temp
  for (int i = 0; i < numberParticles; ++i)
    for(int j = 0; j < 3; ++j)
      particle[i].v[j] *= factor;
}

void applyBC(CVector3& pos)
{
  for (int i = 0; i < 3; ++i)
    pos[i] -= lrint(pos[i] / systemSize[i]) * systemSize[i];
}
double calcVirial(vector<CParticle> &particle, vector<int> &NL, int NLpos[])
{
  double virial = 0;

  //loop through each pair of particles
  for(int i = 0; i < particle.size(); ++i)
    {
      int max_length = NL.size();
      if(i + 1 < particle.size())
        max_length = NLpos[i + 1];
      for(int j = NLpos[i]; j < max_length; ++j)
	{
	  CVector3 distance = particle[i].r - particle[NL[j]].r; //vector between particles
	  applyBC(distance);
	  //calculate virial component
	  if(distance.length() <= 6 * radius)
	  {
	    CVector3 force = calcForce(distance);
	    virial += force.dotProd(distance);
	  }

	}
    }

  return virial;
}
void calcAllForces(vector<CParticle>& particles, vector<int>& NL, int NLpos[])
{
  for(int i = 0; i < particles.size(); ++i)
    {
      int max_length = NL.size();
      if(i + 1 < particles.size())
        max_length = NLpos[i + 1];
      for(int j = NLpos[i]; j < max_length; ++j)
	{
      //cout << i << " " << j << endl;
	  CVector3 distance = particles[i].r - particles[NL[j]].r;
	  applyBC(distance);
	  if(distance.length() <= 6 * radius) //truncate the lennard jones potential at 3 diameters
	    {
	      CVector3 force = calcForce(distance); //calculate the force on the particle
	      //change particle accelerations
	      particles[i].a += force / mass;
	      particles[NL[j]].a -= force / mass;
	    }
	}
    }
}
void zeroMomentum(vector<CParticle> &particle)
{
  double sum[3] = {0, 0, 0};
  for (int i = 0; i < particle.size(); ++i) //loop through particles and calculate sum of x and y velocites
    for(int j = 0; j < 3; ++j)
      sum[j] += particle[i].v[j];
  for (int i = 0; i < particle.size(); ++i)
    for(int j = 0; j < 3; ++j)
      particle[i].v[j] -= sum[j]/particle.size(); //reduce velocity by 1/number of particles of the total sum

}

void calcNeighbourList(vector<CParticle> &particle, vector<int> &NL, int NLpos[])
{
  int count = 0;
  NL.clear(); //clear any old neighbour list
  for(int i = 0; i < particle.size(); ++i)
    {
      NLpos[i] = count; // set position of i's neighbours in neighbour list
      for(int j = i + 1; j < particle.size(); ++j)
	{
	  CVector3 distance = particle[i].r - particle[j].r;
	  applyBC(distance);
	  if(distance.length() <= 6.6 * radius) //if particles are close together
	    {
	      ++count; //increment counter
	      NL.push_back(j); //add particle j to particle i's neighbour list
	    }
	}
    }
}
