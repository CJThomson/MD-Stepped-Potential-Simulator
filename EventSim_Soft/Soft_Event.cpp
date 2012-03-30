//----Program Includes----
#include "Declares.h" //all includes and function declares for event sim
//Physical Properties:
const double density = 0.85;
const double temperature = 1.67; //temperature of the system
double ptail = 32.0 * M_PI * density * density / 9 * (1 / pow(2.3 , 6) - 1.5) / pow(2.3,3);
//Simulation:
int numberParticles = 256; //number of particles
const int numberEvents = 50;
//const double length = pow(numberParticles/density, 1.0/3.0);
const double length = 6.0;
const CVector3 systemSize(length,length,length); //size of the system
double t = 0;
const bool initFile = true; //use an init file instead of random generated values
const bool overwriteInit = false; //create a new init file
std::vector<Steps> steps; //create a vector to store step propeties
bool thermostat = false; //use a thermostat
const double thermoFreq = 0.0001; //update frequency of thermostat
const int thermoOff = 10000;
const int noCells = 3;

//Reduced Unit Definitions
const double mass = 1; //mass of a particlen
const double radius = 0.5; //radius of a particle (set for diameter = 1)
const double sigma = 1;

//Logging:
const int psteps = 50; //frequency of output to file
const int writeOutLog = 1; //level of outLog, 0 = nothing, 1 = event discriptions, 2 = full

//Measuring Properties:
const int noReadings = 200000; //number of readings to take
const int noBins = 1000; //number of radial bins
const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
bool takeMeasurement = false;
double gVal[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over time

//----Time Averages----
double TA_gVal[noBins];
double TA_v = 0;
double TA_U = 0;
double TA_T = 0;
Logger logger; //create an instance of the logger class

std::map<std::pair<int, int>, int> collStep;
using namespace std;

int main()
{
  //variable declarations
  vector<CParticle> particles; //create a vector to store particle info
  vector<vector<eventTimes> > particleEL;
  vector<eventTimes> masterEL; //create a vector to store collision times
  vector<set<int> > neighbourList; //create a vector of sets to hold particles in  neighbour cell
  vector<set<int> > neighbourCell; //create a vector of sets to hold the cells neighbouring each cell
  int readingsTaken = 0;
  cout << "Initialising Simulation with ";
  //Initialise the simulation
  if(initFile)
    initFromFile(particles);
  else
    initialise(particles); // initialise the system
  cout << particles.size() << " particles..." << endl;
  numberParticles = particles.size();

  cout << "Initialising neighbour lists" << endl;
  int cells3 = pow(noCells, 3);
  neighbourList.resize(cells3);
  neighbourCell.resize(cells3);
  generateNeighbourCells(neighbourCell);
  generateNeighbourList(neighbourList, particles);

  cout << "Initialising Steps..." << endl;
  initSteps(); //step up system steps
  if(length * 0.5 < steps[steps.size()-1].step_radius)
    cout << "Warning cell length less than cut-off radius" << endl;
  cout << "Initialising Output Logs..." << endl;
  logger.initialise(Logger::LOCATIONS); //initialise location logger
  logger.initialise(Logger::OUTPUTLOG);

  cout << "Initialising Random Number Generators" << endl;
  //Initialise random number generator
  boost::mt19937 eng;
  //initialise normally distributed random numbers
  boost::normal_distribution <double> nd(0.0, sqrt(3.0 * temperature));
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<double> > var_nor(eng,nd);
  //initialise uniformly distributed random numbers
  boost::uniform_real<double> ud(0.0, 1.0);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > var_uni(eng, ud);

  cout << "Generating event list..." << endl;
  particleEL.resize(particles.size());
  masterEL.resize(particles.size());
  generateEvents(particles, particleEL, masterEL, neighbourList, neighbourCell); //generate the event list
  
  logger.write_Location(particles, 0, systemSize); //write initial values to the log

  cout << "Starting Simulation ..." << endl;
  for(int n = 0; n < numberEvents; ++n)
    {
      bool validEvent = true;
      eventTimes next_event = *min_element(masterEL.begin(), masterEL.end());

      double dt = next_event.collisionTime - t; //find time to next collision

      if((numberEvents - n) <=  noReadings) //take readings this step?
	takeMeasurement = true;

      if(n == thermoOff) //turn thermostat off?
	{
	  thermostat = false;
	  zeroMomentum(particles);
	}

      switch (next_event._type)
	{
	case eventTimes::IP_IN:
	case eventTimes::IP_OUT:
	  {
	    if(next_event.p2coll != particles[next_event.particle2].collNo) //check if collision is valid
	      {
		validEvent = false;
		int particleNo = 0;
		if(next_event.collisionTime == (particleEL[next_event.particle1].front().collisionTime))
		  particleNo = next_event.particle1;
		else
		  particleNo = next_event.particle2;
		updateEvents(particleNo, -1, particles, particleEL, masterEL, neighbourList, neighbourCell);
		logger.outLog << "Invalid Collision between particles " << next_event.particle1 << " & " << next_event.particle2 << endl;
		--n;
	      }
	    else //if a valid event
	      {
		t += dt; //update system time
		updatePositions(particles,dt);
		if(dt < 0)
		  {
		    //cout << "negative dt" << endl;
		  }
		if(writeOutLog == 1)
		  {
		    CVector3 rij = (particles[next_event.particle1].r-particles[next_event.particle2].r);
		    applyBC(rij);
		    logger.outLog << "Collision no: " << n << " between particles: " << next_event.particle1 << " & " << next_event.particle2
				  << " at time = " << t
				  << " at distance of = " << rij.length()
				  << " UEnergy " << calcPotential()
				  << " KEnergy " << calcKinetic(particles)
				  << endl;
		  }
		calcVelocity(particles, next_event);

		updateEvents(next_event.particle1,next_event.particle2,particles, particleEL, masterEL, neighbourList, neighbourCell);
		++particles[next_event.particle1].collNo;
		++particles[next_event.particle2].collNo;
		//apply Andersen Thermostat
		if(thermostat)
		  for(int i = 0; i < numberParticles; ++i)
		    if(var_uni() < thermoFreq)
		      {
			particles[i].v = var_nor() * particles[i].v.normalise();
			updateEvents(i, -1, particles, particleEL, masterEL, neighbourList, neighbourCell);
			logger.outLog << "Thermostat increase on particle " << i << endl;
		      }
	      }

	    break;
	  }
	case eventTimes::SENTINAL:
	  --n;
	  if(writeOutLog == 1)
	    logger.outLog << "Sentinal event for particle: " << next_event.particle1 << " at time = "<< t<< endl;

	  t += dt; //update system time
	  updatePositions(particles, dt);
	  updateEvents(next_event.particle1, -1, particles, particleEL, masterEL, neighbourList, neighbourCell);
	  break;
	case eventTimes::NEIGHBOURCELL:
	  {
	    t += dt; //update system time
	    updatePositions(particles,dt);
	    CParticle p1 = particles[next_event.particle1];
	    neighbourList[p1.cellNo].erase(next_event.particle1); //erase particle from previous cell
	    int sign = -1;
	    if(p1.v[p1.nextCell] >0)
	      sign = 1;
	    int newCell = lrint((p1.cellNo + sign * 1 ) / noCells) * noCells;
	    neighbourList[newCell].insert(p1.particleNo);
	    if(writeOutLog == 1)
	      logger.outLog << "Particle "<< next_event.particle1 << " changed cell from " << p1.cellNo <<
		" to " << newCell <<  " at time = "<< t<< endl;
	    p1.cellNo = newCell;
	    break;
	  }
	case eventTimes::WALL:
	  if(writeOutLog == 1)
	    logger.outLog << "Wall event for particle: " << next_event.particle1 << endl;
	  break;
	case eventTimes::NONE:
	  cout << "Ran out of events!" << endl;
	  exit(1);
	  break;
	default:
	  cout << "Unknown event type! Aborting... " << endl;
	  exit(1);
	  break;
	}
      if((numberEvents - n) ==  noReadings)
	{
	  for(int i = 0; i < particles.size(); ++i)
	    particles[i].r0 = particles[i].r;
	}

      if(n%psteps==0)
	logger.write_Location(particles, t, systemSize);

      if(takeMeasurement && validEvent)
	{
	  ++readingsTaken;
	  if(n%100 == 0)
	    coDiff.push_back(Diffusion(calcDiff(particles, t),t));
	  TA_T += calcTemp(particles);
	  TA_U += calcPotential();
	  if(n%10 == 0)
	    {
	      calcRadDist(particles);
	      for(int i = 0; i < noBins; ++i)
		TA_gVal[i] += gVal[i];
	    }

	}
      if(numberEvents > 100 && validEvent)
	{
	  if(n%(numberEvents / 100) == 0)
	    {
	      cout << n << " of " << numberEvents << " events simulated" <<
		" T: " << calcTemp(particles) <<
		" TE: " << calcPotential() + calcKinetic(particles) << endl;
	    }
	}

    }

  cout << "Simulation Complete" << endl;
  double deltaR = maxR / noBins; //width of each shell

  for(int i = 0; i < noBins; ++i)
    {
      double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
      TA_gVal[i] /= (0.5 * numberParticles * (readingsTaken/10) * volShell * density);
    }
  double freeTime = t * numberParticles / (2 * numberEvents);
  cout << "Time Averages" << endl;
  cout << "Mean free time: " << freeTime <<
    " U: " << TA_U / (readingsTaken * numberParticles) / 2 <<
    " T: " << TA_T / readingsTaken  <<
    " P(short): " << (density) * (TA_T + mass * sigma / (freeTime * 6.0) * TA_v) / readingsTaken <<
    " P(w LR): " << ((density) * (TA_T + mass * sigma / (freeTime * 6.0) * TA_v) / readingsTaken + ptail) << endl;
  cout << "Writing Results..." << endl;
  logger.write_RadDist(TA_gVal,noBins, deltaR); //write radial distribution file
  logger.write_Diff(coDiff); //write coefficient of diffusion file
  logger.write_Location(particles, t, systemSize); //write final values to the log
  if(overwriteInit)
    logger.write_Init(particles); //write final values to log to allow
  cout << "All Tasks Complete" << endl;
  return 0;
}

void applyBC(CVector3& pos)
{
  for (size_t i(0); i < 3; ++i)
    pos[i] -= lrint(pos[i] / systemSize[i])*systemSize[i];
}

int calcCell(CVector3 r)
{
  int cell[3] = {0,0,0};
  double cellSize = length / noCells;
  applyBC(r);
  for(int i = 0; i < 3; ++i)
    cell[i] = floor((r[i] + 0.5 * systemSize[i]) / cellSize);

  return cell[0] * noCells * noCells + cell[1] * noCells + cell[2];
}

double calcCellLeave(CParticle& particle)
{
  double t_min[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  double boundary[3] = {0,0,0};
  double cellLength = length / noCells;
  int cells2 = noCells * noCells;
  int cell[3];
  cell[0] = floor(particle.cellNo / cells2);
  cell[1] = floor((particle.cellNo % cells2) / noCells);
  cell[2] = particle.cellNo % noCells;

  for(int i = 0; i < 3; ++i)
    {
      if(particle.v[i] < 0)
	boundary[i] = cell[i] * cellLength - 0.5 * length;
      else if(particle.v[i] > 0)
	boundary[i] = (cell[i] + 1) * cellLength - 0.5 * length;

      if(particle.v[i] !=0)
	{
	  double distance = boundary[i] - particle.r[i];
	  t_min[i] = distance / fabs(particle.v[i]);
	}
    }
  double tMin =  min(t_min[0], min(t_min[1], t_min[2]));
  if(t_min[0] == tMin)
    particle.nextCell = 0;
  else if (t_min[1] == tMin)
    particle.nextCell = 1;
  else if (t_min[2] == tMin)
    particle.nextCell = 2;
}

double calcCollisionTime(CParticle& particle1, CParticle& particle2, eventTimes::EventType& eventType)
{
  double t_min_out = HUGE_VAL; //set minimum time to infinity
  double t_min_in = HUGE_VAL;
  map<pair<int, int>, int>::iterator it_map;

  int p1 = min(particle1.particleNo, particle2.particleNo);
  int p2 = max(particle1.particleNo, particle2.particleNo);
  CVector3 r12 = particle1.r - particle2.r;
  CVector3 v12 = particle1.v - particle2.v;
  applyBC(r12);
  double r12sqr = r12.dotProd(r12);
  double v12sqr = v12.dotProd(v12);
  double vdotr = v12.dotProd(r12);

  //Output Logging
  if(writeOutLog == 2 )
    {
      logger.outLog << " = = = = NEW COLLISION = = = = " << endl;
      logger.outLog << "Particles involved: " << particle1.particleNo << " & " << particle2.particleNo << endl;
      logger.outLog << "Particle 1 :- r=(" << particle1.r.x << "," <<
	particle1.r.y << "," << particle1.r.z << ") - v=(" <<
	particle1.v.x << "," << particle1.v.y << "," << particle1.v.z << ")" << endl;
      logger.outLog << "Particle 2 :- r=(" << particle2.r.x << "," <<
	particle2.r.y << "," << particle2.r.z << ") - v=(" <<
	particle2.v.x << "," << particle2.v.y << "," << particle2.v.z << ")" << endl;
      logger.outLog << "r12=(" << r12.x << "," <<
	r12.y << "," << r12.z << ") - v12=(" <<
	v12.x << "," << v12.y << "," << v12.z << ")" << endl;
      logger.outLog << "v.r = " << vdotr << endl;
    }

  it_map = collStep.find(pair<int, int> (p1, p2));
  if(it_map == collStep.end()) //if particle is outside steps
    {
      double c = r12sqr - steps[steps.size() - 1].step_radius * steps[steps.size() - 1].step_radius;
      double arg = vdotr * vdotr - v12sqr * c;
      if(arg >= 0) //if particles come near enough to each other
	t_min_in = c / (-vdotr + sqrt(arg));
      if(writeOutLog == 2)
	{
	  logger.outLog << "Particle is outside steps"<< endl;
	  logger.outLog << "c: " << c << " arg: " << arg << endl;
	  logger.outLog << "t_min_in: " << t_min_in << endl;
	}
      eventType = eventTimes::IP_IN;
      if(t_min_in > 0)
	return t + t_min_in;
      else
	return HUGE_VAL;
    }
  else
    {
      if(it_map->second != -1) //if there is an innerstep to interact with
	{
	  if(vdotr < 0)
	    {
	      //calculate collision time inwards
	      double c = r12sqr - steps[it_map->second - 1].step_radius * steps[it_map->second - 1].step_radius;
	      double arg = vdotr * vdotr - v12sqr * c;
	      if(arg >= 0) //if particles come near enough to each other
		t_min_in = c / (-vdotr + sqrt(arg));
	      if(writeOutLog == 2)
		{
		  logger.outLog << "Inward Collision" << endl;
		  logger.outLog << "Step Number: " << it_map->second << endl;
		  logger.outLog << "c: " << c << " arg: " << arg << endl;
		  logger.outLog << "t_min_in: " << t_min_in << endl;
		}
	    }
	}

      //check for any outward steps
      double c = r12sqr - steps[it_map->second].step_radius * steps[it_map->second].step_radius;
      double arg = vdotr * vdotr - v12sqr * c;

      if(arg >= 0) //if particles come near enough to each other
	t_min_out = (sqrt(arg) - vdotr) / v12sqr;

      if(writeOutLog == 2)
	{
	  logger.outLog << "Outward Collision" << endl;
	  logger.outLog << "Step Number: " << it_map->second << endl;
	  logger.outLog << "c: " << c << " arg: " << arg << endl;
	  logger.outLog << "t_min_out: " << t_min_out << endl;
	}
      if(t_min_out < 0)
	{
	  logger.outLog << "INVALID COLLISION TIME" << endl;
	  logger.outLog << "Particles involved: " << particle1.particleNo << " & " << particle2.particleNo << endl;
	  logger.outLog << "Particle 1 :- r=(" << particle1.r.x << "," <<
	    particle1.r.y << "," << particle1.r.z << ") - v=(" <<
	    particle1.v.x << "," << particle1.v.y << "," << particle1.v.z << ")" << endl;
	  logger.outLog << "Particle 2 :- r=(" << particle2.r.x << "," <<
	    particle2.r.y << "," << particle2.r.z << ") - v=(" <<
	    particle2.v.x << "," << particle2.v.y << "," << particle2.v.z << ")" << endl;
	  logger.outLog << "r12=(" << r12.x << "," <<
	    r12.y << "," << r12.z << ") - v12=(" <<
	    v12.x << "," << v12.y << "," << v12.z << ")" << endl;
	  logger.outLog << "Outward Collision" << endl;
	  logger.outLog << "Step Number: " << it_map->second << endl;
	  logger.outLog << "c: " << c << " arg: " << arg << endl;
	  logger.outLog << "v.r = " << vdotr << endl;
	  logger.outLog << "t_min_out: " << t_min_out << endl;
	}

      if(t_min_out < t_min_in || t_min_in < 0)
	{
	  eventType = eventTimes::IP_OUT;
	  return t + t_min_out;
	}
      else
	{
	  eventType = eventTimes::IP_IN;
	  return t + t_min_in;
	}
    }
}

double calcDiff(vector<CParticle>& particle, double time)
{
  double sumDiff = 0;
  for(int i = 0; i < particle.size(); ++i)
    {
      CVector3 distTravelled =  particle[i].r - particle[i].r0;
      sumDiff += distTravelled.dotProd(distTravelled);
    }
  //sumDiff /= (6 * time * particle.size());
  return sumDiff / particle.size();
}

double calcKinetic(vector<CParticle>& particle)
{
  double kinetic = 0;
  for(int i = 0; i < particle.size(); ++i)
    kinetic += particle[i].kineticEnergy();

  return kinetic;
}

double calcPressure(double Temp,double t)
{
  return density * Temp * (1 + Temp / ( 3 * numberParticles * t) * TA_v);
}

double calcPotential()
{
  double potential = 0;
  map<pair<int, int>, int>::iterator i;
  for(i = collStep.begin(); i != collStep.end(); ++i)
    {
      potential += steps[i->second].step_energy;
    }

  return potential;
}

void calcRadDist(vector<CParticle> &particle)
{
  for(int i = 0; i < noBins; ++i)
    gVal[i] = 0; //rezero array

  for(int i = 0; i < particle.size(); ++i)
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

void calcStep(CParticle& particle1, CParticle& particle2)
{
  int p1 = min(particle1.particleNo, particle2.particleNo);
  int p2 = max(particle1.particleNo, particle2.particleNo);
  CVector3 r12 = particle1.r - particle2.r;
  CVector3 v12 = particle1.v - particle2.v;
  applyBC(r12);
  double distance = r12.length();
  for(int i = 1; i < steps.size(); ++i)
    {
      if(distance <= steps[i].step_radius && distance > steps[i - 1].step_radius)
	collStep.insert(pair<pair<int, int>, int> (pair<int, int> (p1, p2), i));
    }
}

double calcSentinalTime(CParticle& particle)
{
  double t_min = HUGE_VAL;
  for(size_t dim = 0; dim < 3; ++dim)
    {
      double vel = fabs(particle.v[dim]);
      double t_sent = 0.25 * (systemSize[dim] - particle.radius) / vel;

      if (vel != 0)
	t_min = min(t_min, t_sent);
    }
  return t + t_min;
}

double calcTemp(vector<CParticle> &particle)
{
  double sum = 0;
  for(unsigned int i=0; i < particle.size(); ++i)
    sum += particle[i].v.dotProd(particle[i].v);
  return mass/(3*particle.size())*sum;
}

void calcVelocity(vector<CParticle>& particle, eventTimes &event)
{
  //variable definitions
  map<pair<int, int>, int>::iterator it_map; //iterator for accessing collision state map
  //set p1 to smaller particle number, and p2 to larger particle number
  int p1 = min(particle[event.particle1].particleNo, particle[event.particle2].particleNo);
  int p2 = max(particle[event.particle1].particleNo, particle[event.particle2].particleNo);

  CVector3 r12 = particle[p1].r - particle[p2].r;
  CVector3 v12 = particle[p1].v - particle[p2].v;
  //apply periodic boundary conditions
  applyBC(r12);
  double vdotr = v12.dotProd(r12.normalise());
  //  =  =  Output Logging  =  =
  if(writeOutLog == 2)
    {
      logger.outLog << " = = = = PROCESSING COLLISION = = = = " << endl;
      logger.outLog << "Particles involved: " << particle[event.particle1].particleNo << " & " << particle[event.particle2].particleNo << endl;
      logger.outLog << "Particle 1 :- r=(" << particle[event.particle1].r.x << "," <<
	particle[event.particle1].r.y << "," << particle[event.particle1].r.z << ") - v=(" <<
	particle[event.particle1].v.x << "," << particle[event.particle1].v.y << "," << particle[event.particle1].v.z << ")" << endl;
      logger.outLog << "Particle 2 :- r=(" << particle[event.particle2].r.x << "," <<
	particle[event.particle2].r.y << "," << particle[event.particle2].r.z << ") - v=(" <<
	particle[event.particle2].v.x << "," << particle[event.particle2].v.y << "," << particle[event.particle2].v.z << ")" << endl;
      logger.outLog << "r12=(" << r12.x << "," <<
	r12.y << "," << r12.z << ") - v12=(" <<
	v12.x << "," << v12.y << "," << v12.z << ")" << endl;
    }
  switch(event._type)
    {
    case eventTimes::IP_IN:
      {
	it_map = collStep.find(pair<int,int>(p1, p2)); //find collision state of particles
	double dU = 0;
	if(it_map == collStep.end() ) //if no collision state found then particles must be outside outer step
	  {
	    dU = -steps[steps.size() - 1].step_energy; //energy is the outermost step height
	    collStep.insert(pair<pair<int, int>, int> (pair<int, int> (p1, p2), steps.size())); //insert pair into collStep map
	    it_map = collStep.find(pair<int,int>(p1, p2)); //set iterator to new value
	  }
	else //if not outside then energy change is the difference in step heights
	  dU = steps[it_map->second].step_energy - steps[it_map->second - 1].step_energy;

	if(writeOutLog == 2) //if writing full outlog
	  {
	    logger.outLog << "IP_IN" << endl;
	    logger.outLog << "dU: " << dU << endl;
	    logger.outLog << "v.r: " << vdotr << endl;
	    if(it_map != collStep.end())
	      logger.outLog << "Step being processed: " << it_map->second << endl;
	  }

	if((vdotr * vdotr + 4.0 / mass * dU) > 0) //if particles go over the step
	  {
	    double A = -0.5 / mass * (vdotr + sqrt(vdotr * vdotr +  4.0 / mass * dU)); //change in momentum
	    if(writeOutLog == 2) {logger.outLog << "A: " << A << endl;}
	    if(takeMeasurement) //if taking measurements
	      {
		TA_v += A / mass;
	      }

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();
	    if(it_map->second != 0)
	      --(it_map->second); //move particles in one step

	  }
	else //if bounce occurs
	  {
	    if(takeMeasurement) //if taking measurements
	      TA_v -= vdotr;
	    particle[p1].v -= vdotr * r12.normalise();
	    particle[p2].v += vdotr * r12.normalise();
	  }
	break;
      }
    case eventTimes::IP_OUT:
      {

	it_map = collStep.find(pair<int, int> (p1, p2));
	double dU = 0;
	if(it_map->second != steps.size() - 1)
	  dU = steps[it_map->second].step_energy - steps[it_map->second + 1].step_energy; //step particle is going to - step particle is on
	else
	  dU = steps[it_map->second].step_energy;
	if(writeOutLog ==2)
	  {
	    logger.outLog << "IP_OUT" << endl;
	    logger.outLog << "Step being processed: " << it_map->second << endl;
	    logger.outLog << "dU: " << dU << endl;
	    logger.outLog << "v.r: " << vdotr << endl;
	  }
	if((vdotr * vdotr + 4.0 / mass * dU) > 0) //if particles go over the step
	  {
	    double A = -0.5 / mass * (vdotr - sqrt(vdotr * vdotr +  4.0 / mass * dU)); //change in momentum
	    if(writeOutLog == 2) {logger.outLog << "A: " << A << endl;}
	    if(takeMeasurement) //if taking measurements
	      {
		TA_v += A / mass;
	      }

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();
	    ++(it_map->second); //move particles in one step
	    if(it_map->second == steps.size()) //if particle moves outside final step delete it's entry
	      collStep.erase(it_map);
	  }
	else //if bounce occurs
	  {
	    if(takeMeasurement) //if taking measurements
	      TA_v -= vdotr;
	    particle[p1].v -= vdotr * r12.normalise();
	    particle[p2].v += vdotr * r12.normalise();
	  }
	break;
      }
    }
  logger.outLog << endl;
}

void correctVelocity(vector<CParticle> &particle)
{
  zeroMomentum(particle);
  double temp = calcTemp(particle); //calculate temp of system
  double factor = sqrt(temperature / temp); //calculated correction factor
  for(int i = 0; i < numberParticles; ++i)
    for(int j = 0; j < 3; ++j)
      particle[i].v[j] *= factor; //scale velocity
}

void initialise(vector<CParticle> &particle)
{
  int n = ceil(pow(numberParticles/4,(double)1/3)); //find cubic root of number of particles
  double scale = 1;
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
      double v2= double(rand()%200)/100-1;
      double s2= u*u + v*v;
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
      switch (j)
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
		  particle.push_back(CParticle(init_r, init_v, radius, mass, particle.size()));
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

void initSteps()
{
  //Steps from Chepela et al, Case 6
  steps.push_back(Steps(0.80,66.74)); //step 0
  steps.push_back(Steps(0.85,27.55)); //step 1
  steps.push_back(Steps(0.90, 10.95)); //step 3
  steps.push_back(Steps(0.95, 3.81)); //step 4
  steps.push_back(Steps(1.00, 0.76)); //step 5
  steps.push_back(Steps(1.05, -0.47)); //step 6
  steps.push_back(Steps(1.25,-0.98)); //step 7
  steps.push_back(Steps(1.45,-0.55)); //step 8
  steps.push_back(Steps(1.75,-0.22)); //step 9
  steps.push_back(Steps(2.30,-0.06)); //step 10
}

void generateEvents(vector<CParticle>& particle, //vector of particles
		    vector<vector<eventTimes> >& pEvents,
		    vector<eventTimes>& events, //master event list
		    vector<set<int> >& NL, //the particles in each cell
		    vector<set<int> >& NC) //the cells neighbouring each cell
{
  set<int>::iterator it_NL;
  set<int>::iterator it_NC;
  for(size_t i = 0; i < particle.size(); ++i)
    {
      double t_min = calcSentinalTime(particle[i]);
      pEvents[i].push_back(eventTimes(t_min, i, -1, 0, eventTimes::SENTINAL));
      //t_min = calcCellLeave(particle[i]);
      //pEvents[i].push_back(eventTimes(t_min,i ,-1, eventTimes::NEIGHBOURCELL));

      /*for(it_NC = NC[particle[i].cellNo].begin(); it_NC != NC[particle[i].cellNo].end(); ++it_NC)
	for(it_NL = NL[*it_NC].begin(); it_NL != NL[*it_NC].end(); ++it_NL)*/
      for(int j = 0; j < particle.size();++j)
	{
	  calcStep(particle[i], particle[j]);
	  eventTimes::EventType eventType;
	  double t_min = calcCollisionTime(particle[i], particle[j], eventType);
	  pEvents[i].push_back(eventTimes(t_min, i, j,  particle[j].collNo, eventType));
	}
      make_heap(pEvents[i].begin(), pEvents[i].end(), greater<eventTimes>());
      events[i] = pEvents[i].front();
    }
}

void generateNeighbourCells(vector<set<int> >& NC)
{
  int cells2 = noCells * noCells;

  for(int i = 0; i < NC.size(); ++i)
    {
      int xmain = floor(i / cells2);
      int ymain = floor((i % cells2) / noCells);
      int zmain = i % noCells;
      for(int dx = -1; dx < 2; ++dx)
	for(int dy = -1; dy < 2; ++dy)
	  for(int dz = -1; dz < 2; ++dz)
	    {
	      int xnew = xmain + dx;
	      if(xnew < 0) {xnew += noCells;}
	      else if(xnew >= noCells) {xnew -= noCells;}

	      int ynew = ymain + dy;
	      if(ynew < 0) {ynew += noCells;}
	      else if(ynew >= noCells) {ynew -= noCells;}

	      int znew = zmain + dz;
	      if(znew < 0) {znew += noCells;}
	      else if(znew >= noCells) {znew -= noCells;}

	      NC[i].insert(xnew * cells2 + ynew * noCells + znew);
	    }
    }
}

void generateNeighbourList(vector<set<int> >& NL, vector<CParticle>& particle)
{
  for(int i = 0; i < numberParticles; ++i)
    {
      int cell = calcCell(particle[i].r);
      NL[cell].insert(particle[i].particleNo);
      particle[i].cellNo = cell;
    }
}

void updateEvents(int particle1, int particle2,
		  vector<CParticle>& particle,
		  //vector<priority_queue<eventTimes, vector<eventTimes>, greater<eventTimes> > >& pEvents, //each particle eventlist
		  vector<vector<eventTimes> >& pEvents,
		  vector<eventTimes> &events, //master event list
		  vector<set<int> >& NL, //the particles in each cell
		  vector<set<int> >& NC) //the cells neighbouring each cell
{
  set<int>::iterator it_NL;
  set<int>::iterator it_NC;
  pEvents[particle1].clear();

  double t_min = calcSentinalTime(particle[particle1]);
  pEvents[particle1].push_back(eventTimes(t_min, particle1, -1, 0, eventTimes::SENTINAL));
  //t_min = calcCellLeave(particle[particle1]);
  //pEvents[particle1].push_back(eventTimes(t_min, particle1 ,-1, eventTimes::NEIGHBOURCELL));
  /*for(it_NC = NC[particle[particle1].cellNo].begin(); it_NC != NC[particle[particle1].cellNo].end(); ++it_NC)
    for(it_NL = NL[*it_NC].begin(); it_NL !=NL[*it_NC].end(); ++it_NL)*/
  for(int j = 0; j < particle.size();++j)
    {
      eventTimes::EventType eventType;
      double t_min = calcCollisionTime(particle[particle1], particle[j], eventType);
      pEvents[particle1].push_back(eventTimes(t_min, particle1, j, particle[j].collNo ,eventType));
    }
  make_heap(pEvents[particle1].begin(), pEvents[particle1].end(), greater<eventTimes>());
  events[particle1] = pEvents[particle1].front();

  if(particle2 != -1)
    {
      pEvents[particle2].clear();

      double t_min = calcSentinalTime(particle[particle2]);
      pEvents[particle2].push_back(eventTimes(t_min, particle2, -1, 0, eventTimes::SENTINAL));
      //t_min = calcCellLeave(particle[particle2]);
      //pEvents[particle2].push_back(eventTimes(t_min, particle2 ,-1, eventTimes::NEIGHBOURCELL));

      /*for(it_NC = NC[particle[particle2].cellNo].begin(); it_NC != NC[particle[particle2].cellNo].end(); ++it_NC)
	for(it_NL = NL[*it_NC].begin(); it_NL !=NL[*it_NC].end(); ++it_NL)*/
      for(int j = 0; j < particle.size();++j)
	{
	  eventTimes::EventType eventType;
	  double t_min = calcCollisionTime(particle[particle2], particle[j], eventType);
	  pEvents[particle2].push_back(eventTimes(t_min, particle2, j, particle[j].collNo, eventType));
	}
      make_heap(pEvents[particle2].begin(), pEvents[particle2].end(), greater<eventTimes>());
      events[particle2] = pEvents[particle2].front();
    }
}

void updatePositions(vector<CParticle>& particle, double dt)
{
  for(int i = 0; i < particle.size(); ++i)
    particle[i].r += particle[i].v * dt;
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
