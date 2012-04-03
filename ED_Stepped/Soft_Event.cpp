//----Program Includes----
#include "Declares.h" //all includes and function declares for event sim
//Physical Properties:
const double density = 0.85;
const double temperature = 1.67; //temperature of the system
double ptail = 32.0 * M_PI * density * density / 9 * (1 / pow(2.3 , 6) - 1.5) / pow(2.3,3);
//Simulation:
int numberParticles = 256; //number of particles
const int numberEvents = 50000;
int eventCount = 0;
const double length = pow(numberParticles/density, 1.0 / 3.0);
//const double length = 6.0;
const CVector3 systemSize(length,length,length); //size of the system
double t = 0;
const bool initFile = false; //use an init file instead of random generated values
const bool overwriteInit = false; //create a new init file
std::vector<Steps> steps; //create a vector to store step propeties
bool thermostat = false; //use a thermostat
const double thermoFreq = 0.0001; //update frequency of thermostat
const int thermoOff = 10000;
const int noCells = 3;

//Reduced Unit Definitions
const double mass = 1; //mass of a particlen
const double radius = 0.5; //radius of a particle (set for diameter = 1)

//Logging:
const int psteps = 50; //frequency of output to file
const int writeOutLog = 1;//level of outLog, 0 = nothing, 1 = event discriptions, 2 = full

//Measuring Properties:
const int noReadings = 200000; //number of readings to take
const int noBins = 1000; //number of radial bins
const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
bool takeMeasurement = false;
double gVal[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over

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
  srand(10);
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
    cout << "Warning system size less than cut-off radius" << endl;
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
  //for(vector<CParticle>::iterator p1 = particles.begin(); p1 != particles.end(); ++p1)
  for(size_t i = 0; i < particles.size(); ++i)
    getEvent(particles[i], particles, particleEL, masterEL, neighbourList, neighbourCell);
  
  logger.write_Location(particles, 0, systemSize); //write initial values to the log

  cout << "Starting Simulation ..." << endl;
  for(;eventCount < numberEvents;)
    {
      bool validEvent = true;
      eventTimes next_event = *min_element(masterEL.begin(), masterEL.end());

      double dt = next_event.collisionTime - t; //find time to next collision

      if((numberEvents - eventCount) <=  noReadings) //take readings this step?
	takeMeasurement = true;

      if(eventCount == thermoOff) //turn thermostat off?
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
		//cerr << "Invalid Collision between particles :" << next_event.particle1 << " & " << next_event.particle2 << endl;
		validEvent = false;
		getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
		getEvent(particles[next_event.particle2], particles, particleEL, masterEL, neighbourList, neighbourCell);
		logger.outLog << "Invalid Collision between particles: " << next_event.particle1 << " & " << next_event.particle2 << endl;
	      }
	    else //if a valid event
	      {
		t += dt; //update system time
		if(writeOutLog >= 1)
		  {
		    updatePosition(particles[next_event.particle1]);
		    updatePosition(particles[next_event.particle2]);
		    CVector3 rij = (particles[next_event.particle1].r-particles[next_event.particle2].r);
		    applyBC(rij);
		    //cerr << "Collision between particles: " << next_event.particle1 << " & " << next_event.particle2 << endl;
		    logger.outLog << "Collision no: " << eventCount << " between particles: " << next_event.particle1 << " & " << next_event.particle2
				  << " at time = " << t
				  << " at distance of = " << rij.length()
				  << " UEnergy " << calcPotential()
				  << " KEnergy " << calcKinetic(particles)
				  << endl;
		  }
		calcVelocity(particles, next_event);		
		++particles[next_event.particle1].collNo;
		++particles[next_event.particle2].collNo;

		getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
		getEvent(particles[next_event.particle2], particles, particleEL, masterEL, neighbourList, neighbourCell);
		++eventCount;
	      }

	    break;
	  }
	case eventTimes::SENTINAL:
	  if(writeOutLog >= 1)
	    logger.outLog << "Sentinal event for particle: " << next_event.particle1 << " at time = "<< t<< endl;

	  t += dt; //update system time
	  getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
	  break;
	case eventTimes::NEIGHBOURCELL:
	  {
	    t += dt; //update system time
	    CParticle p1 = particles[next_event.particle1];
	    neighbourList[p1.cellNo].erase(next_event.particle1); //erase particle from previous cell
	    int sign = -1;
	    if(p1.v[p1.nextCell] >0)
	      sign = 1;
	    int newCell = lrint((p1.cellNo + sign * 1 ) / noCells) * noCells;
	    neighbourList[newCell].insert(p1.particleNo);
	    if(writeOutLog >= 1)
	      logger.outLog << "Particle "<< next_event.particle1 << " changed cell from " << p1.cellNo
			    << " to " << newCell <<  " at time = "<< t<< endl;
	    p1.cellNo = newCell;
	    break;
	  }
	case eventTimes::WALL:
	  if(writeOutLog >= 1)
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
      if((numberEvents - eventCount) ==  noReadings)
	{
	  for(int i = 0; i < particles.size(); ++i)
	    particles[i].r0 = particles[i].r;
	}

      if(eventCount%psteps==0)
	logger.write_Location(particles, t, systemSize);

      if(takeMeasurement && validEvent)
	{
	  ++readingsTaken;

	  if(eventCount%100 == 0)
	    coDiff.push_back(Diffusion(calcDiff(particles, t),t));
	  TA_T += calcTemp(particles);
	  TA_U += calcPotential();
	  if(eventCount%10 == 0)
	    {
	      calcRadDist(particles);
	      for(int i = 0; i < noBins; ++i)
		TA_gVal[i] += gVal[i];
	    }

	}
#if(numberEvents >= 100)
      if(n%(numberEvents / 100) == 0 && validEvent)
	{
	  cout << n << " of " << numberEvents << " events simulated"
	       << " T: " << calcTemp(particles)
	       << " TE: " << calcPotential() + calcKinetic(particles) << endl;
	}
#endif
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
  cout << "Mean free time: " << freeTime
       << " U: " << TA_U / (readingsTaken * numberParticles) / 2
       << " T: " << TA_T / readingsTaken
       << " P: " << (density) * (TA_T + mass / (freeTime * 6.0) * TA_v) / readingsTaken + ptail
       << endl;
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
  for(size_t i(0); i < 3; ++i)
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

  for(size_t i(0); i < 3; ++i)
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

  return tMin;
}

double calcCollisionTime(CParticle& particle1, CParticle& particle2, eventTimes::EventType& eventType)
{
  double t_min_out = HUGE_VAL; //set minimum time to infinity
  double t_min_in = HUGE_VAL;
  map<pair<int, int>, int>::const_iterator it_map;

  CVector3 r12 = particle1.r - particle2.r;
  CVector3 v12 = particle1.v - particle2.v;
  applyBC(r12);
  double r12sqr = r12.dotProd(r12);
  double v12sqr = v12.dotProd(v12);
  double vdotr = v12.dotProd(r12);

  //Output Logging
  if(writeOutLog >= 2 )
    {
      logger.outLog << " = = = = NEW COLLISION = = = = " << endl;
      logger.outLog << "Particles involved: " << particle1.particleNo << " & " << particle2.particleNo << endl;
      logger.outLog << "Particle 1 :- r=(" << particle1.r.x << ","
		    << particle1.r.y << "," << particle1.r.z << ") - v=("
		    << particle1.v.x << "," << particle1.v.y << "," << particle1.v.z << ")" << endl;
      logger.outLog << "Particle 2 :- r=(" << particle2.r.x << ","
		    << particle2.r.y << "," << particle2.r.z << ") - v=("
		    << particle2.v.x << "," << particle2.v.y << "," << particle2.v.z << ")" << endl;
      logger.outLog << "r12=(" << r12.x << ","
		    << r12.y << "," << r12.z << ") - v12=("
		    << v12.x << "," << v12.y << "," << v12.z << ")" << endl;
      logger.outLog << "v.r = " << vdotr << endl;
    }

  it_map = collStep.find(pair<int, int> (min(particle1.particleNo, particle2.particleNo), 
					 max(particle1.particleNo, particle2.particleNo)));
  if (it_map == collStep.end()) //if particle is outside steps
    {
      double c = r12sqr - steps[steps.size() - 1].step_radius * steps[steps.size() - 1].step_radius;
      double arg = vdotr * vdotr - v12sqr * c;
      if((vdotr < 0) && (arg >= 0)) //if particles come near enough to each other
	t_min_in = c / (-vdotr + sqrt(arg));
      if(writeOutLog >= 2)
	{
	  logger.outLog << "Particle is outside steps"<< endl;
	  logger.outLog << "c: " << c << " arg: " << arg << endl;
	  logger.outLog << "t_min_in: " << t_min_in << endl;
	}
      eventType = eventTimes::IP_IN;
      return (t_min_in > 0) ? t + t_min_in : HUGE_VAL;
    }
  else
    {
      if(it_map->second != 0) //if there is an innerstep to interact with
	{
	  if(vdotr < 0)
	    {
	      //calculate collision time inwards
	      double c = r12sqr - steps[it_map->second - 1].step_radius * steps[it_map->second - 1].step_radius;
	      double arg = vdotr * vdotr - v12sqr * c;
	      if(arg >= 0) //if particles come near enough to each other
		t_min_in = c / (-vdotr + sqrt(arg));
	      if(writeOutLog >= 2)
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
      else //There's been a numerical error, so just collide at their nearest point
	t_min_out = - vdotr / v12sqr;

      if(writeOutLog >= 2)
	{
	  logger.outLog << "Outward Collision" << endl;
	  logger.outLog << "Step Number: " << it_map->second << endl;
	  logger.outLog << "c: " << c << " arg: " << arg << endl;
	  logger.outLog << "t_min_out: " << t_min_out << endl;
	}
      if(t_min_out < 0)
	{
	  //cerr << "Negative dt between particles: " << p1 << " and " << p2 << endl;
	  logger.outLog << "INVALID COLLISION TIME (event=" << eventCount << ")" << endl;
	  logger.outLog << "Particles involved: " << particle1.particleNo << " & " << particle2.particleNo << endl;
	  logger.outLog << "Particle 1 :- r=(" << particle1.r.x << ","
			<< particle1.r.y << "," << particle1.r.z << ") - v=("
			<< particle1.v.x << "," << particle1.v.y << "," << particle1.v.z << ")" << endl;
	  logger.outLog << "Particle 2 :- r=(" << particle2.r.x << ","
			<< particle2.r.y << "," << particle2.r.z << ") - v=("
			<< particle2.v.x << "," << particle2.v.y << "," << particle2.v.z << ")" << endl;
	  logger.outLog << "r12=(" << r12.x << ","
			<< r12.y << "," << r12.z << ") - v12=("
			<< v12.x << "," << v12.y << "," << v12.z << ")" << endl;
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

double calcDiff(vector<CParticle>& particles, double time)
{
  double sumDiff = 0;
  for(vector<CParticle>::iterator particle = particles.begin(); particle != particles.end(); ++particle)
    {
      CVector3 distTravelled =  particle->r - particle->r0;
      sumDiff += distTravelled.dotProd(distTravelled);
    }
  return sumDiff / particles.size();
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
    potential += steps[i->second].step_energy;

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
	if(distance.length() < maxR)
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
  for(int i = 0; i < steps.size(); ++i)
    {
      if(distance <= steps[i].step_radius)
	{
	  collStep.insert(pair<pair<int, int>, int> (pair<int, int> (p1, p2), i));
	  break;
	}
    }
}

double calcSentinalTime(CParticle& particle)
{
  double t_min = HUGE_VAL;
  for(size_t dim (0); dim < 3; ++dim)
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

void calcVelocity(vector<CParticle>& particle, eventTimes& event)
{
  //variable definitions
  map<pair<int, int>, int>::iterator it_map; //iterator for accessing collision state map

  updatePosition(particle[event.particle1]);
  updatePosition(particle[event.particle2]);


  //set p1 to smaller particle number, and p2 to larger particle number
  int p1 = min(event.particle1, event.particle2);
  int p2 = max(event.particle1, event.particle2);

  CVector3 r12 = particle[p1].r - particle[p2].r;
  CVector3 v12 = particle[p1].v - particle[p2].v;
  //apply periodic boundary conditions
  applyBC(r12);
  double vdotr = v12.dotProd(r12.normalise());
  //  =  =  Output Logging  =  =
  if(writeOutLog >= 2)
    {
      logger.outLog << " = = = = PROCESSING COLLISION = = = = " << endl;
      logger.outLog << "Particles involved: " << particle[event.particle1].particleNo
		    << " & " << particle[event.particle2].particleNo << endl;
      logger.outLog << "Particle 1 :- r=(" << particle[event.particle1].r.x << "," 
		    << particle[event.particle1].r.y << "," << particle[event.particle1].r.z << ") - v=("
		    << particle[event.particle1].v.x << "," << particle[event.particle1].v.y
		    << "," << particle[event.particle1].v.z << ")" << endl;
      logger.outLog << "Particle 2 :- r=(" << particle[event.particle2].r.x << ","
		    << particle[event.particle2].r.y << "," << particle[event.particle2].r.z
		    << ") - v=(" << particle[event.particle2].v.x << ","
		    << particle[event.particle2].v.y << "," << particle[event.particle2].v.z << ")"
		    << endl;
      logger.outLog << "r12=(" << r12.x << "," << r12.y << "," << r12.z << ") - v12=("
		    << v12.x << "," << v12.y << "," << v12.z << ")" << endl;
    }
  switch(event._type)
    {
    case eventTimes::IP_IN:
      {
	it_map = collStep.find(pair<int,int>(p1, p2)); //find collision state of particles
	double dU = 0;
	if(it_map == collStep.end()) //if no collision state found then particles must be outside outer step
	  dU = -steps.back().step_energy; //energy is the outermost step height
	else //if not outside then energy change is the difference in step heights
	  dU = steps[it_map->second].step_energy - steps[it_map->second - 1].step_energy;

	if(writeOutLog >= 2) //if writing full outlog
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
	    if(writeOutLog >= 2) 
	      {
		logger.outLog << "A: " << A << endl;
		logger.outLog << "Well Capture" << endl;
	      }
	    if(takeMeasurement) //if taking measurements
	      TA_v += r12.length() * A / mass;

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();

	    if(it_map == collStep.end()) //if no collision state found then particles must be outside outer step
	      collStep.insert(pair<pair<int, int>, int>(pair<int, int>(p1, p2), steps.size() - 1)); //insert pair into collStep map
	    else
	      --(it_map->second); //move particles in one step
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >= 2)
	      logger.outLog << "Well Bounce" << endl;
	    if(takeMeasurement) //if taking measurements
	      TA_v -= r12.length() * vdotr;
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

	if(writeOutLog >=2)
	  {
	    logger.outLog << "IP_OUT" << endl;
	    logger.outLog << "Step being processed: " << it_map->second << endl;
	    logger.outLog << "dU: " << dU << endl;
	    logger.outLog << "v.r: " << vdotr << endl;
	  }
	if((vdotr * vdotr + 4.0 / mass * dU) > 0) //if particles go over the step
	  {
	    double A = -0.5 / mass * (vdotr - sqrt(vdotr * vdotr +  4.0 / mass * dU)); //change in momentum
	    if(writeOutLog >= 2) 
	      {
		logger.outLog << "A: " << A << endl;
		logger.outLog << "Well Release" << endl;
	      }
	    if(takeMeasurement) //if taking measurements
	      TA_v += r12.length() * A / mass;

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();
	    if(it_map->second == steps.size() -1) //if particle is leaving outermost step
	      collStep.erase(it_map);
	    else
	      ++(it_map->second); //move particles in one step
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >=2)
	      logger.outLog << "Well Bounce" << endl;
	    if(takeMeasurement) //if taking measurements
	      TA_v -= r12.length() * vdotr;
	    particle[p1].v -= vdotr * r12.normalise();
	    particle[p2].v += vdotr * r12.normalise();
	  }
	break;
      }
    default:
      std::cerr << "Unhandled calcVelocities type" << std::endl;
      std::exit(1);
    }

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

      double R = sqrt(-1 * log(s) / s);
      double R2 = sqrt(-1 * log(s2) / s2);
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
  //steps.push_back(Steps(1.00,0.8)); //step 0
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

void generateNeighbourList(vector<set<int> >& NL, vector<CParticle>& particles)
{
  for(vector<CParticle>::iterator particle = particles.begin(); particle != particles.end(); ++particle)
    {
      int cell = calcCell(particle->r);
      NL[cell].insert(particle->particleNo);
      particle->cellNo = cell;
    }
}

void getEvent(CParticle& p1,
	      vector<CParticle>& particles,
	      vector<vector<eventTimes> >& pEvents,
	      vector<eventTimes> &events, //master event list
	      vector<set<int> >& NL, //the particles in each cell
	      vector<set<int> >& NC) //the cells neighbouring each cell )
{
  set<int>::iterator it_NL;
  set<int>::iterator it_NC;
  pEvents[p1.particleNo].clear();
  updatePosition(p1);
  double t_min_sent = calcSentinalTime(p1); 
  pEvents[p1.particleNo].push_back(eventTimes(t_min_sent, p1.particleNo, -1, 0, eventTimes::SENTINAL));
  //t_min = calcCellLeave(particle[particle1]);
  //pEvents[particle1].push_back(eventTimes(t_min, particle1 ,-1, eventTimes::NEIGHBOURCELL));
  /*for(it_NC = NC[particle[particle1].cellNo].begin(); it_NC != NC[particle[particle1].cellNo].end(); ++it_NC)
    for(it_NL = NL[*it_NC].begin(); it_NL !=NL[*it_NC].end(); ++it_NL)*/
  //for(vector<CParticle>::iterator p2 = particles.begin(); p2 != particles.end(); ++p2)
  for(size_t j = 0; j < particles.size(); ++j)
    {
      if (p1.particleNo == j) continue;
      /*if(j > p1.particleNo)
	cerr << " p1 " << p1.particleNo << " p2 " << j << endl;*/
      updatePosition(particles[j]);
      eventTimes::EventType eventType;
      double t_min_coll = calcCollisionTime(p1, particles[j], eventType);
      pEvents[p1.particleNo].push_back(eventTimes(t_min_coll, p1.particleNo, j, particles[j].collNo, eventType));
      /*if (p1.particleNo == p2->particleNo) break;
      if(j > p1.particleNo)
	cerr << "yay p1 " << p1.particleNo << " p2 " << p2->particleNo << endl;
      eventTimes::EventType eventType;
      double t_min_coll = calcCollisionTime(p1, *p2, eventType);
      pEvents[p1.particleNo].push_back(eventTimes(t_min_coll, p1.particleNo, p2->particleNo, p2->collNo, eventType));
      if (p1.particleNo == 111 || p1.particleNo == 62)
	if(p2->particleNo == 111 || p2->particleNo == 62)
	cerr << "t_min_coll" << t_min_coll << endl;*/
    }
  events[p1.particleNo] = *min_element(pEvents[p1.particleNo].begin(), pEvents[p1.particleNo].end());
}
void updatePosition(CParticle& particle)
{
  double delta_t = t - particle.updateTime;
  if(delta_t != 0)
    {
      particle.r += delta_t * particle.v;
      particle.updateTime = t; 
    }
}

void zeroMomentum(vector<CParticle> &particles)
{
  double sum[3] = {0, 0, 0};
  for (vector<CParticle>::iterator particle = particles.begin(); particle < particles.end(); ++particle) //loop through particles and calculate sum of x and y velocites
    for(size_t j (0); j < 3; ++j)
      sum[j] += particle->v[j];
  
  for (vector<CParticle>::iterator particle = particles.begin(); particle < particles.end(); ++particle) //loop through particles and calculate sum of x and y velocites
    for(size_t j (0); j < 3; ++j)
      particle->v[j] -= sum[j]/particles.size(); //reduce velocity by 1/number of particles of the total sum

}
