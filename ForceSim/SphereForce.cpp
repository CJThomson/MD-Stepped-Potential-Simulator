//----Program Includes----
#include "Declares.h" //declarations for this program

//Physical Properties:
double density = 0.85;
double temperature = 4.6; //temperature of the system
//Simulation:
int numberParticles = 864; //number of particles
const int simTime = 50000; //length of the Simulation
const double dt = 0.005; //length of ticme interval
double length = pow(numberParticles / density, 1.0 / 3.0);
const double r_cut = 3.0;
CVector3 systemSize(length, length, length); //size of the system
const bool initFile = false; //use an init file
const bool overwriteInit = false; //create a new initfile
const int NL_update = 10;
int number_of_runs = 10;
//const double ptail = -16 * M_PI * density * density / (3.0 * pow(r_cut, 3)) * (1.0 - 2.0 / (3.0 * pow(r_cut,6))); 
//const double utail = -8 * M_PI * density / (3.0 * pow(r_cut, 3)) * (1.0 - 1.0 / (3.0 * pow(r_cut,6))); 

const double ptail = 0;
const double utail = 0;
//Thermostat
bool thermostat = true; //use a thermostat
const double thermoFreq = 0.05; //update frequency of thermostat
const int thermoOff = simTime;

//Reduced Unit Definitions:
const double mass = 1.0; //mass of a particle
const double radius = 0.5; //radius of a particle (set for diameter = 1)
const double epsilon = 1.0; //minimum energy of Lennard Jones Potential
const double sigma = 1.0; //distance for Lennard Jones root

//Logging:
Logger logger; //create instance of the logger class
const int out_interval = 20; //frequency of output to file
const int diff_interval = simTime;
const int rdf_interval = 10;
const int sample_interval = 10;
const bool writeLoc = false;

//Measuring Properties:
const int startSampling = 20000; //number of readings to take
double readingTime = 0;
const int noBins = 600; //number of radial bins
const double maxR = 3.0;
//const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
bool suppressOutput = true;

double gVal[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over time

//----Time Averages----
double TA_gVal[noBins];
double TA_Temp = 0;
double TA_U = 0;
double TA_Virial = 0;
double TA_Temp2 = 0;
double TA_U2 = 0;
double TA_Virial2 = 0;

using namespace std;
int main()
{
  string input;
  cout << " - - - Force-Driven Simulator - - - " << endl;
  while(true)
    {
      cout << "input: ";
      cin >> input;
      if(input == "exit")
	exit(0);
      else if(input == "start")
	{
	  vector<Results> results;
	  results.resize(number_of_runs);
	  cout << endl;
	  for(size_t i (0); i < number_of_runs; ++i)
	    {
	      cout << "\rRun number: " << i + 1 << " of " << number_of_runs << endl;
	      resetSim();
	      runSimulation(results, i);
	    }
	  Results avgResults;
	  for(vector<Results>::iterator result = results.begin(); result != results.end(); ++result)
	    avgResults += *result;
	  avgResults *= 1.0 / number_of_runs;
	  cout << "Average Results" << endl;
	  cout << "Temperature: " << avgResults.temperature << endl;
	  cout << "Pressure: " << avgResults.pressure << endl;
	  cout << "Potential Energy: " << avgResults.potential << endl;
	  logger.write_Results(results, density, temperature);
	  //After Simulation
	  double deltaR = maxR / noBins; //width of each shell
	  for(size_t i(0); i < noBins; ++i)
	    {
	      double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
	      TA_gVal[i]/=(0.5 * numberParticles * (simTime - startSampling) * number_of_runs / rdf_interval * volShell * density);
	    }
	  logger.write_RadDist(TA_gVal,noBins, deltaR, density, temperature); //write radial gdistribution file
	}
      else if(input == "temperature")
	cin >> temperature;
      else if(input == "density")
	cin >> density;
      else if(input == "particles")
	cin >> numberParticles;
      else if(input == "runs")
	cin >> number_of_runs;
      else
	cout << "Invalid input" << endl;
    }
}
void resetSim()
{
  length = pow(numberParticles / density, 1.0 / 3.0);
  for(size_t i(0); i < 3; ++i)
    systemSize[i] = length; 
  thermostat = true;
  readingTime = 0;
  coDiff.clear(); 
  TA_Temp = 0;
  TA_U = 0;
  TA_Virial = 0;
  TA_Temp2 = 0;
  TA_U2 = 0;
  TA_Virial2 = 0;

}

void runSimulation(vector<Results>& results, size_t runNumber)
{
  //variable declarations

  vector<CParticle> particles;
  int noReadings = 0;

  cout << "Initialising random number generator...";
  CRandom RNG;
  RNG.seed();
  //Initialise the simulation
  if(initFile)
    initFromFile(particles); //initialise the system from file
  else
    initialise(particles, RNG); // initialise the system
  numberParticles = particles.size();

  cout << "\rInitialising particles" << flush;



  cout << "\rInitialising Log files..." << flush;
  if(writeLoc) { logger.initialise(Logger::LOCATIONS); } //initialise location logger
  if(writeLoc) { logger.write_Location(particles, 0, systemSize); } //write initial values to log

  cout << "\rInitialising neighbour lists..." << flush;
  vector<int> neighbourList;
  int listPos[particles.size()];
  calcNeighbourList(particles, neighbourList, listPos);
  calcAllForces(particles, neighbourList, listPos); //calculate all forces on system


  cout << "\rStarting Simulation with " << numberParticles
       << " at a density of " << density << " and a temperature of " << temperature << " in a box of length " << length << endl;
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
	logger.write_Location(particles, t, systemSize); //write locations to output file

      if(timeStep % NL_update == 0)
	calcNeighbourList(particles, neighbourList, listPos);

      if(timeStep % 100 == 0) //write a screen output every 100 timesteps
	{
	  double pshort = (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3)) + ptail;
      	  cout << "\r" << fixed << setprecision(2) << (double) timeStep / simTime * 100 << "% complete" << "\t"
	       << " t: " << setprecision(2) << t << "\t"
	       << " T: " << setprecision(4) << calcTemp(particles) << "\t"
	       << " P: " << setprecision(3) << pshort << "\t"
	       << " U: " << setprecision(3) << calcPotential(particles,neighbourList, listPos) << "\t"
	       << " TE: " << calcPotential(particles, neighbourList, listPos) + calcKinetic(particles)
	       << flush;
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

      if(timeStep > startSampling)
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
	      double temp_virial = calcVirial(particles, neighbourList, listPos);
	      double temp_kinetic = calcKinetic(particles);
	      TA_Virial += temp_virial;
	      TA_Virial2 += temp_virial * temp_virial;
	      double temp_T = calcTemp(particles);
	      TA_Temp += temp_T; //add current temp to TA
	      TA_Temp2 += temp_T * temp_T;
	      double temp_potential = calcPotential(particles, neighbourList, listPos);
	      TA_U += temp_potential; //add current potential energy to TA
	      TA_U2 += temp_potential * temp_potential;
	    }
	}

    }

  cout << "\rSimulation Complete" << flush;
  //After Simulation
   cout << "\rwrite output files..." << flush;
  logger.write_Location(particles, simTime, systemSize);
  logger.write_Diff(coDiff); //write coefficient of diffusion file
  if(overwriteInit)
    logger.write_Init(particles);

  double E_Temp = TA_Temp * sample_interval / noReadings;
  double E_Press = density * E_Temp
    + TA_Virial * density  * sample_interval/ (3.0 * numberParticles * noReadings);
  double E_U = TA_U * sample_interval / noReadings;
  results[runNumber] = Results(E_Temp, E_Press + ptail, E_U +utail);
  //output Time Averaged system properties
  if(!suppressOutput)
    {
      cout << setprecision(5) << "Time Averages:" << endl;
      // = Tempature
      double E_Temp = TA_Temp * sample_interval / noReadings;
      double temp_sd = sqrt(fabs((TA_Temp2 * sample_interval / noReadings) - E_Temp * E_Temp));
      cout << "Temp: " << E_Temp << "(" << temp_sd << ")" << endl;
      
      // = Pressure
      double pshort = (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3));
      double E_Press = density * E_Temp
	+ TA_Virial * density  * sample_interval/ (3.0 * numberParticles * noReadings);
      double press_sd = sqrt( pow( density * temp_sd, 2) 
			      + pow( density * sample_interval / (3.0 * numberParticles * noReadings), 2));
      cout << "P: " << E_Press << "(" << press_sd << ")" << endl;
      cout << "P+LR: " << E_Press + ptail << endl;
      
      // = Potential Energy
      double E_U = TA_U * sample_interval / noReadings;
      double U_sd = sqrt(fabs(TA_U2 * sample_interval / noReadings - E_U * E_U));
      cout << "U: " << E_U << "(" << U_sd << ")" << endl;
      cout << "U+LR: " << E_U + utail << endl;
    }
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
  for(size_t i (0); i < numberParticles; ++i)
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

double calcDiff(vector<CParticle>& particles, double time)
{
  double sumDiff = 0;
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    {
      CVector3 distTravelled = particle->r - particle->r0;
      sumDiff += distTravelled.dotProd();
    }

  return sumDiff / particles.size();
}

double calcTemp(vector<CParticle> &particles)
{
  double sum = 0;
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    sum += particle->v.dotProd();
  return mass / ( 3 * particles.size()) * sum;
}

void calcRadDist(vector<CParticle> &particles)
{
  for(size_t i (0); i < noBins; ++i)
    gVal[i] = 0; //rezero array
  
  for(it_part p1 = particles.begin(); p1 != particles.end(); ++p1)
    for(it_part p2 = p1 + 1; p2 != particles.end(); ++p2)
      {
	CVector3 distance = p1->r - p2->r;
	applyBC(distance);
	if(distance.length() < maxR)
	  {
	    int index = floor(distance.length() * noBins / maxR);
	    ++gVal[index];
	  }
      }
}

CVector3 calcForce(CVector3 rij)
{
  return 24.0 * epsilon / sigma * (2.0 * pow((sigma / rij.length()), 13) - pow((sigma / rij.length()), 7)) * rij.normalise();
}
double calcKinetic(vector<CParticle> &particles)
{
  double kinetic = 0;
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    kinetic += particle->kineticEnergy();
  return kinetic;
}

double calcPotential(vector<CParticle> &particles, vector<int> &NL, int NLpos[])
{
  double potential = 0;
  //loop through each pair of particles
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    {
      int max_length = NL.size();
      if(particle->particleNo + 1 < particles.size())
        max_length = NLpos[particle->particleNo + 1];
      for(size_t j = NLpos[particle->particleNo]; j < max_length; ++j)
	{
	  CVector3 distance = particle->r - particles[NL[j]].r; //vector between particles
	  applyBC(distance);
	  if(distance.length() <= 6 * radius)
	    potential += 4.0 * epsilon * (pow(sigma / distance.length(), 12) - pow(sigma / distance.length(), 6));
	}
    }
  return potential / particles.size();
}

void correctVelocity(vector<CParticle> &particles)
{
  zeroMomentum(particles); //zero linear momentum of the system
  double temp = calcTemp(particles); //calculate temperature of system
  double factor = sqrt(temperature/temp); //calculate correction factor to set system to correct temp
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    for(size_t dim (0); dim < 3; ++dim)
      particle->v[dim] *= factor;
}

void applyBC(CVector3& pos)
{
  for (size_t dim (0); dim < 3; ++dim)
    pos[dim] -= lrint(pos[dim] / systemSize[dim]) * systemSize[dim];
}
double calcVirial(vector<CParticle> &particles, vector<int> &NL, int NLpos[])
{
  double virial = 0;

  //loop through each pair of particles
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    {
      int max_length = NL.size();
      if(particle->particleNo + 1 < particles.size())
        max_length = NLpos[particle->particleNo + 1];
      for(size_t j = NLpos[particle->particleNo]; j < max_length; ++j)
	{
	  CVector3 distance = particle->r - particles[NL[j]].r; //vector between particles
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

  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    {
      int max_length = NL.size();
      if(particle->particleNo + 1 < particles.size())
        max_length = NLpos[particle->particleNo + 1];
      for(size_t j = NLpos[particle->particleNo]; j < max_length; ++j)
	{
	  CVector3 distance = particle->r - particles[NL[j]].r; //vector between particles
	  applyBC(distance);
	  if(distance.length() <= 6 * radius)
	    {
	      CVector3 force = calcForce(distance); //calculate the force on the particle
	      //change particle accelerations
	      particle->a += force / mass;
	      particles[NL[j]].a -= force / mass;
	    }
	}
    }
}

void zeroMomentum(vector<CParticle> &particles)
{
  double sum[3] = {0, 0, 0};
  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    for(size_t dim = 0; dim < 3; ++dim)
      sum[dim] += particle->v[dim];

  for(it_part particle = particles.begin(); particle != particles.end(); ++particle)
    for(size_t dim = 0; dim < 3; ++dim)
      particle->v[dim] -= sum[dim]/particles.size(); //reduce velocity by 1/number of particles of the total sum

}

void calcNeighbourList(vector<CParticle> &particles, vector<int> &NL, int NLpos[])
{
  int count = 0;
  NL.clear(); //clear any old neighbour list
  for(it_part p1 = particles.begin(); p1 != particles.end(); ++p1)
    {
      NLpos[p1->particleNo] = count; // set position of i's neighbours in neighbour list
      for(it_part p2 = p1 + 1; p2 != particles.end(); ++p2)
	{
	  CVector3 distance = p1->r - p2->r;
	  applyBC(distance);
	  if(distance.length() <= 6.6 * radius) //if particles are close together
	    {
	      ++count; //increment counter
	      NL.push_back(p2->particleNo); //add particle j to particle i's neighbour list
	    }
	}
    }
}
