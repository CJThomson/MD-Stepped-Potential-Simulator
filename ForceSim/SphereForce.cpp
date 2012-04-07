//----Program Includes----
#include "Declares.h" //declarations for this program

//Physical Properties:
const double density = 0.85;
const double temperature = 1.3; //temperature of the system

//Simulation:
const int numberParticles =256; //number of particles
const double simTime = 550; //length of the Simulation
const double dt = 0.005; //length of time interval
const double length = pow(numberParticles / density, 1.0 / 3.0);
const CVector3 systemSize(length, length, length); //size of the system
const bool initFile = false; //use an init file
const bool overwriteInit = false; //create a new initfile
bool thermostat = true; //use a thermostat
const double thermoFreq = 0.05; //update frequency of thermostat
//Reduced Unit Definitions:
const double mass = 1; //mass of a particle
const double radius = 0.5; //radius of a particle (set for diameter = 1)
const double epsilon = 1; //minimum energy of Lennard Jones Potential
const double sigma = 1; //distance for Lennard Jones root

//Logging:
const int psteps = 20; //frequency of output to file
const bool writeLoc = false;

//Measuring Properties:
const int noReadings = 30000; //number of readings to take
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
  //Initialise the simulation
  if(initFile)
    initFromFile(particles); //initialise the system from file
  else
    initialise(particles); // initialise the system

  vector<int> neighbourList;
  int listPos[particles.size()];

  if(writeLoc) { log.initialise(Logger::LOCATIONS); } //initialise location logger
  //log.initialise(Logger::BOLTZMANN_H); //initialise logging of Boltzmann H value
  calcNeighbourList(particles, neighbourList, listPos);
  calcAllForces(particles, neighbourList, listPos); //calculate all forces on system
  if(writeLoc) { log.write_Location(particles, 0, systemSize); } //write initial values to log

  //Initialise random number generator
  boost::mt19937 eng;
  //initialise normally distributed random numbers
  boost::normal_distribution <double> nd(0.0, sqrt(3 * temperature));
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > var_nor(eng,nd);
  //initialise uniformly distributed random numbers
  boost::uniform_real<double> ud(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > var_uni(eng, ud);

  for (double t = dt; t < simTime; t += dt)
    {
      //move all particles forward
      for(int i = 0; i < numberParticles; ++i)
	{
	  //move particles in space using Velocity Verlet
	  particles[i].r += dt * (particles[i].v + particles[i].a * dt * 0.5);
	  particles[i].v += dt * particles[i].a * 0.5; //calculate part of velcity
	  particles[i].a = CVector3(0,0,0); //zero acceleration
	}

      calcAllForces(particles, neighbourList, listPos); //calculate forces (and accelerations)

      for(int i = 0; i < numberParticles; ++i)
	particles[i].v += dt * particles[i].a * 0.5; //calculate remainder of velocity

      //apply Andersen Thermostat
      if(thermostat)
	{
	  for(int i = 0; i < numberParticles; ++i)
	    {
	      if(var_uni() < thermoFreq)
		{
		  /*for(int j = 0; j < 3; ++j)
		    particles[i].v[j] = var_nor();*/
		    particles[i].v = particles[i].v.normalise() * var_nor();
		}
	    }
	}

      int readNo = (t - simTime) / dt + noReadings;
      //update particle
      if(int(t/dt)%psteps==0)
	{
	  log.write_Location(particles, t, systemSize); //write locations to output file
	  //log.write_BoltzH(t, calcBoltzmannH(particles));
	  if(readNo >0)
	    {
	      double D = calcDiff(particles, t); //calculate value of diffusion coefficient
	      coDiff.push_back(Diffusion(D,t)); //add value to vector
	    }

	}
      if(int(t/dt)%5 == 0)
	calcNeighbourList(particles, neighbourList, listPos);
      if(int(t/dt)%100==0) //write a screen output every 100 timesteps
	{
	  double rc=3.0;
	  double irc3 = 1.0/(rc*rc*rc), irc6 = irc3*irc3;
	  double ptail = 32*M_PI*density*density/9*(irc6 - 1.5)*irc3;

	  double pshort = (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3)) + ptail;

      	  cout << fixed << setprecision(2) << t/simTime * 100 << "% complete" << "\t" <<
	    " t: " << setprecision(2) << t << "\t" <<
	    " T: " << setprecision(4) << calcTemp(particles) << "\t" <<
	    " P: " << setprecision(3) << pshort + ptail  <<
	    " U: " << setprecision(3) << calcPotential(particles,neighbourList, listPos)  <<
	    endl;
	}

      if(readNo + 20000 >= 0 && thermostat)
      {
        thermostat = false;  //turn off thermostat just prior to taking readings
        zeroMomentum(particles); //zero linear momentum in the system

      }

      if(readNo == 0)
	{
	  readingTime = t;
	  for(int i = 0; i < particles.size(); ++i)
	    particles[i].r0 = particles[i].r;
	}

      if(readNo > 0 && readNo%4 == 0) //only take a reading during final noReadings timesteps and only 1 every 4 steps
	{
	  calcRadDist(particles); //calculate radial distribution
	  for(int i = 0; i < noBins; ++i)
	    TA_gVal[i]+=gVal[i]; //calculate running total

	  TA_Pressure += (2.0 * calcKinetic(particles) + calcVirial(particles, neighbourList, listPos)) / (3.0 * pow(length,3)); //add current pressure to TA
	  TA_Temp += calcTemp(particles); //add current temp to TA
	  TA_U += calcPotential(particles, neighbourList, listPos); //add current potential energy to TA
	}
    }
  //After Simulation
  double deltaR = maxR / noBins; //width of each shell
  for(int i = 0; i < noBins; ++i)
    {
      double volShell = 4.0 / 3.0 * PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
      TA_gVal[i]/=(0.5 * numberParticles * noReadings * 0.25 * volShell * density);
    }

  log.write_Location(particles, simTime, systemSize);
  log.write_RadDist(TA_gVal,noBins, deltaR); //write radial gdistribution file
  log.write_Diff(coDiff); //write coefficient of diffusion file
  if(overwriteInit)
    log.write_Init(particles);

  double rc=3.0;
  double irc3 = 1.0/(rc*rc*rc), irc6 = irc3*irc3;
  double ptail = 32*M_PI*density*density/9*(irc6 - 1.5)*irc3;

  //output Time Averaged system properties
  cout << "Time Averages:- Temp: " << TA_Temp / (noReadings * 0.25) <<
    " Pressure(Total): " << TA_Pressure / (noReadings * 0.25) + ptail  <<
    " U: " << TA_U / (numberParticles * noReadings * 0.25 ) << endl;
  return 0;
}
void initFromFile (vector<CParticle> &particle)
{
  ifstream initLog ("init86.dat");
  string line;
  int strpos = 1;
  while(initLog.good() )
    {
      int counter = 0; //rezero counter
      CVector3 init_r;
      CVector3 init_v;
      getline(initLog, line); //get next line of log
      vector<CParticle> p;
      for(int i = 0; i < line.size(); ++i)
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

void initialise(vector<CParticle>& particle)
{
  int n = ceil(pow(numberParticles/4,(double)1/3)); //find cubic root of number of particles
  double a = length / (n);
  int j = 0, x = 0, y = 0, z = 0;
  for(int i = 0; i < numberParticles; ++i)
    {
      //create two uniformly distrubuted random numbers
      double u = double(rand()%200)/100-1;
      double v = double(rand()%200)/100-1;
      double s = u*u + v*v;
      while(s>=1) //if u^2 + v^2 is creater than zero, ie point lies outside a circle radius 1
        {
	  //generate some new random numbers
	  u = double(rand()%200)/100-1;
	  v = double(rand()%200)/100-1;
	  s = u*u + v*v;
        }
      double u2 = double(rand()%200)/100-1;
      double v2 = double(rand()%200)/100-1;
      double s2 = u*u + v*v;
      while(s2>=1) //if u^2 + v^2 is creater than zero, ie point lies outside a circle radius 1
        {
	  //generate some new random numbers
	  u2 = double(rand()%200)/100-1;
	  v2 = double(rand()%200)/100-1;
	  s2 = u2*u2 + v2*v2;
        }

      double R = sqrt(-1*log(s)/s);
      double R2 = sqrt(-1*log(s2)/s2);
      CVector3 location;
      switch (j) //set particles in a FCC structure
        {
        case 0:
	  location = CVector3(x * a - length * 0.5, y * a - length * 0.5, z * a - length * 0.5);
	  break;
        case 1:
	  location = CVector3(x* a - length * 0.5, (y + 0.5) * a - length * 0.5, (z + 0.5) * a - length * 0.5);
	  break;
        case 2:
	  location = CVector3((x + 0.5) * a - length * 0.5, y * a - length * 0.5, (z + 0.5) * a - length * 0.5);
	  break;
        case 3:
	  location = CVector3((x + 0.5) * a - length * 0.5, (y + 0.5) * a - length * 0.5, z * a - length * 0.5);
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

      CVector3 velocity(u*R, v*R, u2*R2); //assign a random unit vector to the velocity
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
	  if(distance.length() <= maxR)
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
double calcBoltzmannH(vector<CParticle> &particle)
{
  const double vmax = 20;
  const int noBins = 1000;
  const double dv = vmax / noBins;
  double H[3] = {0, 0, 0};
  int bin[3][noBins];

  for(int i = 0; i < noBins; ++i)
    {
      for(int k = 0; k < 3; ++k)
	bin[k][i] = 0;
      for(int j = 0; j < particle.size(); ++j)
	{
	  double v = i * 2 * dv - vmax;
	  for(int k = 0; k < 3; ++k)
	    {
	      if(particle[j].v[k] > v - dv && particle[j].v[k] < v + dv)
		++bin[k][i];
	    }
	}
      for(int k = 0; k < 3; ++k)
	{
	  if(bin[k][i] !=0)
	    {
	      double f = (double) bin[k][i] / particle.size();
	      H[k] += f * log(f) * dv;
	    }


	}

    }

  return (H[0] + H[1] + H[2]) / 3.0;
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
