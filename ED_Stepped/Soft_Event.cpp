//----Program Includes----
#include "Declares.h" //all includes and function declares for event sim
//Physical Properties:
const double density = 0.85;
const double temperature = 1.34; //temperature of the system
//Simulation:
int numberParticles = 256; //number of particles
const int numberEvents = 5e+6;
int eventCount = 0;
const double length = pow(numberParticles/density, 1.0 / 3.0);
//const double length = 6.0;
const CVector3 systemSize(length,length,length); //size of the system
double t = 0;
const bool initFile = false; //use an init file instead of random generated values
const bool overwriteInit = false; //create a new init file
std::vector<Steps> steps; //create a vectOr to store step propeties
const int noCells = 3;

//Thermostat:
bool thermostat = true; //use a thermostat
const size_t thermoFreq = 200; //frequency thermostat rate is updated
size_t thermoCount = 0; //counter to change thermostat rate
size_t thermoLastUpdate = 0;
double thermoMeanFreeTime = 0.0005;
//const int thermoOff =3.5e+6;
const int thermoOff = numberEvents;
double thermoSetting = 0.05;

//Reduced Unit Definitions
const double mass = 1; //mass of a particle
const double radius = 0.5; //radius of a particle (set for diameter = 1)

//Logging:
const int psteps = 50; //frequency of output to file
const int writeOutLog = 0;//level of outLog, 0 = nothing, 1 = event discriptions, 2 = full

//Measuring Properties
const int startSampling = 4e+6; //number of readings to take
const int sample_interval = 1;
const int rdf_interval = 100;
const int diff_interval = 20;
double lastSampleTime = 0;
double startSampleTime = -1;
const int noBins = 500; //number of radial bins
const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
double gVal[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over

//----Time Averages----
double TA_gVal[noBins];
double TA_v = 0;
double TA_U = 0;
double TA_T = 0;
double TA_p = 0;
double TA_v2 = 0;
double TA_U2 = 0;
double TA_T2 = 0;
double TA_p2 = 0;
double TA_tavg2 = 0;
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
  cout << "Initialising Random Number Generators" << endl;
  CRandom RNG;

  cout << "Initialising ";
  //Initialise the simulation
  if(initFile)
    initFromFile(particles);
  else
    initialise(particles, RNG ); // initialise the system
  cout << particles.size() << " particles..." << endl;
  numberParticles = particles.size();
  for(it_particle p1 = particles.begin(); p1 != particles.end(); ++p1)
    for(it_particle p2 = p1 + 1; p2 != particles.end(); ++p2)
      calcStep(*p1,*p2);

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


  cout << "Generating event list..." << endl;
  particleEL.resize(particles.size()); 
  masterEL.resize(particles.size() + 1);//one event for each particle plus thermostat
  for(it_particle p1 = particles.begin(); p1 != particles.end(); ++p1)
    getEvent(*p1, particles, particleEL, masterEL, neighbourList, neighbourCell);

  int particleNo = -1;
  double t_min_thermo = calcThermoTime(RNG, particleNo);
  if(!thermostat)
    t_min_thermo = HUGE_VAL;
  if(particleNo == -1)
    {
      cerr << "calcThermoTime has not returned a valid particle No" << endl;
      exit(1);
    }
  masterEL[particles.size()] = eventTimes(t_min_thermo, particleNo, -1, -1, eventTimes::THERMOSTAT);
  
  logger.write_Location(particles, 0, systemSize); //write initial values to the log

  cout << "Starting Simulation ..." << endl;
  for(;eventCount < numberEvents;)
    {
      bool collEvent = false;
      double temp_pFlux = 0;
      eventTimes next_event = *min_element(masterEL.begin(), masterEL.end());

      double dt = next_event.collisionTime - t; //find time to next collision

      if(eventCount == thermoOff) //turn thermostat off?
	{
	  thermostat = false;
	  masterEL.back().collisionTime = HUGE_VAL;
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
		getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
		getEvent(particles[next_event.particle2], particles, particleEL, masterEL, neighbourList, neighbourCell);
		logger.outLog << "Invalid Collision between particles: " << next_event.particle1 << " & " << next_event.particle2 << endl;
	      }
	    else //if a valid event
	      {
		collEvent = true;
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
		temp_pFlux = calcVelocity(particles, next_event);		
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
	    int newCell = calcNewCell(p1);
	    neighbourList[newCell].insert(p1.particleNo);
	    if(writeOutLog >= 1)
	      logger.outLog << "Particle "<< next_event.particle1 << " changed cell from " << p1.cellNo
			    << " to " << newCell <<  " at time = "<< t<< endl;
	    p1.cellNo = newCell;
	    getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
	    break;
	  }
	case eventTimes::THERMOSTAT:
	  {
	    if(writeOutLog >= 1)
	      logger.outLog << "Thermostat event for particle: " << next_event.particle1 << " at time = "<< t<< endl;
	    t += dt; //update system time
	    runThermostat(particles[next_event.particle1], RNG, masterEL);
	    ++particles[next_event.particle1].collNo;
	    getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
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

      if(collEvent)
	{
	  if(eventCount % psteps == 0) //output file logging
	    logger.write_Location(particles, t, systemSize);

	  if(eventCount == startSampling)
	    {
	      // store particle positions when starting to take readings
	      for(it_particle particle = particles.begin(); particle != particles.end(); ++particle)
		particle->r0 = particle->r;
	      startSampleTime = t;
	    }

	  if(eventCount > startSampling)
	    {
	      if(startSampleTime <0) {cerr<< "ERROR" << endl; exit(1);}
	      if(eventCount % sample_interval == 0)
		{
		  ++readingsTaken;
		  double deltaT = t - lastSampleTime;
		  TA_tavg2 += deltaT * deltaT;
		  // = Temperature
		  double temp_temperature = calcTemp(particles);
		  TA_T += temp_temperature * deltaT;
		  TA_T2 += temp_temperature * temp_temperature * deltaT;

		  // = Potential
		  double temp_potential = calcPotential();
		  TA_U += temp_potential * deltaT;
		  TA_U2 += temp_potential * temp_potential * deltaT;

		  // = Momentum Flux
		  TA_v += temp_pFlux * deltaT;
		  TA_v2 += temp_pFlux * temp_pFlux * deltaT;

		  // = Pressure
		  double temp_pressure = (density * temp_temperature + temp_pFlux * mass * density / (3.0 * particles.size() * deltaT));
		  TA_p += temp_pressure * deltaT;
		  TA_p2 += temp_pressure * temp_pressure * deltaT;
		  lastSampleTime = t;
		}
	      if(eventCount % diff_interval == 0)
		coDiff.push_back(Diffusion(calcDiff(particles, t),t));

	      if(eventCount % rdf_interval == 0)
		{
		  calcRadDist(particles);
		  for(size_t i = 0; i < noBins; ++i)
		    TA_gVal[i] += gVal[i];
		}
	      
	    }
	  if(eventCount % (int) ceil(numberEvents / 1000) == 0)
	    {
	      cout << "\r" << eventCount << " of " << numberEvents << " events simulated"
		   << " T: " << calcTemp(particles)
		   << " TE: " << calcPotential() + calcKinetic(particles) << flush;
	    }
	}
    }

  cout << endl << "Simulation Complete" << endl;
  double deltaR = maxR / noBins; //width of each shell

  for(int i = 0; i < noBins; ++i)
    {
      double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
      TA_gVal[i] /= (0.5 * numberParticles * (readingsTaken/ rdf_interval) * volShell * density);
    }
  //Output time averages
  cout << "Time Averages:" << endl;
  {
    // = Mean free time
    double freeTime = t * numberParticles / (2 * eventCount);
    double freeTime2 = TA_tavg2 * numberParticles / (2 * readingsTaken);
    cout << "Mean free time: " << freeTime << " (" << freeTime2 - freeTime * freeTime << ")" << endl;
    
    // = Potential Energy
    double E_pot = TA_U / t;
    double E_pot2 = TA_U2 / t;
    cout << " U: " << E_pot << "(" << E_pot2 - E_pot * E_pot << ")" << endl;
    
    // = Temperature
    double E_temp = TA_T / t;
    double E_temp2 = TA_T2 / t;
    cout << " T: " << E_temp << "(" << E_temp2 - E_temp * E_temp << ")" << endl;
    
    // = Pressure
    double E_press = TA_p / t;
    double E_press2 = TA_p2 / t;
    cout << " P: " << E_press << "(" << E_press2 - E_press * E_press << ")" << endl;
    
    // = Momentum flux
    double E_mf = TA_v / t;
    double E_mf2 = TA_v2 / t;
    cout << " <r.v> " << E_mf << "(" << E_mf2 - E_mf * E_mf << ")" << endl << endl;;
  }
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

int calcNewCell(CParticle& particle)
{
  cerr << "here" << endl;
  int cells2 = noCells * noCells;
  int cell[3];
  cell[0] = floor(particle.cellNo / cells2);
  cell[1] = floor((particle.cellNo % cells2) / noCells);
  cell[2] = particle.cellNo % noCells;
  int sign = (particle.v[particle.nextCell] < 0) ? -1 : 1;
  cell[particle.nextCell] += sign;
  cerr << cell[particle.nextCell] << endl;
  cerr << lrint(cell[particle.nextCell] / noCells) << endl;
  cell[particle.nextCell] -= floor((double)cell[particle.nextCell] / noCells) * noCells;
  cerr << cell[particle.nextCell] << endl;
  return cell[0] * noCells * noCells + cell[1] * noCells + cell[2];
}
double calcCellLeave(CParticle& particle)
{
  updatePosition(particle);
  CVector3 location = particle.r;
  applyBC(location);
  double t_min[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
  double boundary[3] = {0,0,0};
  double cellLength = length / noCells;
  int cells2 = noCells * noCells;
  int cell[3];
  cell[0] = floor(particle.cellNo / cells2);
  cell[1] = floor((particle.cellNo % cells2) / noCells);
  cell[2] = particle.cellNo % noCells;
   cerr << " = = = START = = =" << endl;
  cerr << particle.cellNo << endl;
  cerr << cellLength << endl;
  cerr << "(" << cell[0] << ", " << cell[1] << ", " << cell[2] << ")" << endl;
  for(size_t i(0); i < 3; ++i)
    {
      if(particle.v[i] < 0)
	boundary[i] = cell[i] * cellLength - 0.5 * length;
      else if(particle.v[i] > 0)
	boundary[i] = (cell[i] + 1) * cellLength - 0.5 * length;

      if(particle.v[i] !=0)
	{
	  double distance = boundary[i] - location[i];
	  t_min[i] = distance / particle.v[i];
	}
    }
  cerr << "b=(" << boundary[0] << ", " << boundary[1] << ", " << boundary[2] << ")" << endl;
  cerr << "l=(" << location[0] << ", " << location[1] << ", " << location[2] << ")" << endl;
  cerr << "v=(" << particle.v[0] << ", " << particle.v[1] << ", " << particle.v[2] << ")" << endl;
  cerr << "v=(" << t_min[0] << ", " << t_min[1] << ", " << t_min[2] << ")" << endl;
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

double calcKinetic(vector<CParticle>& particles)
{
  double kinetic = 0;
  for(vector<CParticle>::iterator particle = particles.begin(); particle != particles.end(); ++particle)
    kinetic += particle->kineticEnergy();

  return kinetic / particles.size();
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

  return potential / numberParticles;
}

void calcRadDist(vector<CParticle> &particles)
{
  for(size_t i (0); i < noBins; ++i)
    gVal[i] = 0; //rezero array

  for(it_particle p1 = particles.begin(); p1 != particles.end(); ++p1)
    {
      updatePosition(*p1);
      for(it_particle p2 = p1 + 1; p2 != particles.end(); ++p2)
	{
	  updatePosition(*p2);
	  CVector3 distance = p1->r - p2->r;
	  applyBC(distance);
	  if(distance.length() < maxR)
	    {
	      int index = floor(distance.length() * noBins/ maxR);
	      if(index < 0 || index >= noBins)
		{
		  cerr << "ERROR: Invalid index is calcRadDist: " << index << endl;
		  exit(1);
		}
	      ++gVal[index];
	    }
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

double calcTemp(vector<CParticle> &particles)
{
  double sum = 0;
  for(it_particle particle = particles.begin(); particle != particles.end(); ++particle)
    sum += particle->v.dotProd();
  return mass / (3 * particles.size()) * sum;
}

double calcVelocity(vector<CParticle>& particle, eventTimes& event)
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

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();

	    if(it_map == collStep.end()) //if no collision state found then particles must be outside outer step
	      collStep.insert(pair<pair<int, int>, int>(pair<int, int>(p1, p2), steps.size() - 1)); //insert pair into collStep map
	    else
	      --(it_map->second); //move particles in one step

	    return r12.length() * A / mass;
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >= 2)
	      logger.outLog << "Well Bounce" << endl;
	    particle[p1].v -= vdotr * r12.normalise();
	    particle[p2].v += vdotr * r12.normalise();
	    return -vdotr;
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

	    //update particle velocities
	    particle[p1].v += A / mass * r12.normalise();
	    particle[p2].v -= A / mass * r12.normalise();
	    if(it_map->second == steps.size() -1) //if particle is leaving outermost step
	      collStep.erase(it_map);
	    else
	      ++(it_map->second); //move particles in one step
	    return r12.length() * A / mass;
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >=2)
	      logger.outLog << "Well Bounce" << endl;
	    particle[p1].v -= vdotr * r12.normalise();
	    particle[p2].v += vdotr * r12.normalise();
	    return -vdotr;
	  }
	break;
      }
    default:
      std::cerr << "Unhandled calcVelocities type" << std::endl;
      std::exit(1);
    }

}

void correctVelocity(vector<CParticle> &particles)
{
  zeroMomentum(particles);
  double temp = calcTemp(particles); //calculate temp of system
  double factor = sqrt(temperature / temp); //calculated correction factor
  for(it_particle particle = particles.begin(); particle != particles.end(); ++particle)
    for(size_t dim (0); dim < 3; ++dim)
      particle->v[dim] *= factor; //scale velocity
}

void initialise(vector<CParticle> &particles, CRandom& RNG)
{
  int n = ceil(pow(numberParticles/4,(double)1/3)); //find cubic root of number of particles
  double a = length / n;
  int j = 0, x = 0, y = 0, z = 0;
  for(size_t i (0); i < numberParticles; ++i)
    {
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
      if(j == 0)
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

      CVector3 velocity; //assign a random unit vector to the velocity
      for(size_t dim(0); dim < 3; ++dim)
	velocity[dim] = RNG.var_normal();

      particles.push_back(CParticle(location, velocity, radius, mass, i));
    }

  correctVelocity(particles);
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
  //steps.push_back(Steps(1.00,500)); //step 0
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

  for(size_t i (0); i < NC.size(); ++i)
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
  for(it_particle particle = particles.begin(); particle != particles.end(); ++particle)
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
  set<int>::iterator p2;
  set<int>::iterator it_NC;
  pEvents[p1.particleNo].clear();
  updatePosition(p1);
  double t_min_sent = calcSentinalTime(p1); 
  pEvents[p1.particleNo].push_back(eventTimes(t_min_sent, p1.particleNo, -1, 0, eventTimes::SENTINAL));
  /*uble t_min_NL = calcCellLeave(p1);
  pEvents[p1.particleNo].push_back(eventTimes(t_min_NL, p1.particleNo ,-1,0, eventTimes::NEIGHBOURCELL));
  for(it_NC = NC[p1.cellNo].begin(); it_NC != NC[p1.cellNo].end(); ++it_NC)
  for(p2 = NL[*it_NC].begin(); p2 !=NL[*it_NC].end(); ++p2)*/
  for(it_particle p2 = particles.begin(); p2 != particles.end(); ++p2)
    {
      if (p1.particleNo == p2->particleNo) continue;
      updatePosition(*p2);
      eventTimes::EventType eventType;
      double t_min_coll = calcCollisionTime(p1, *p2, eventType);
      pEvents[p1.particleNo].push_back(eventTimes(t_min_coll, p1.particleNo, p2->particleNo, p2->collNo, eventType));
      /*if (p1.particleNo == *p2) continue;
      updatePosition(particles[*p2]);
      eventTimes::EventType eventType;
      double t_min_coll = calcCollisionTime(p1, particles[*p2], eventType);
      pEvents[p1.particleNo].push_back(eventTimes(t_min_coll, p1.particleNo, *p2, particles[*p2].collNo, eventType));
      */
    }
  events[p1.particleNo] = *min_element(pEvents[p1.particleNo].begin(), pEvents[p1.particleNo].end());
}
double calcThermoTime(CRandom& RNG, int& particleNo)
{
  particleNo = RNG.var_uniformInt(0, numberParticles - 1); //generate particle to 'collide' with thermostat
  //cerr << particleNo << endl;
  double t_min = -thermoMeanFreeTime * log(RNG.var_01());
  //cerr << t_min << endl;
  return t + t_min;

}

void runThermostat(CParticle& particle, CRandom& RNG, vector<eventTimes> &events)
{
  CVector3 oldv = particle.v;
  ++thermoCount;
  if(thermoCount > thermoFreq)
    {
      if((eventCount - thermoLastUpdate != 0)) //check divisior is not zero
	{
	  thermoMeanFreeTime *= (double) thermoCount 
	    / ((eventCount - thermoLastUpdate) * thermoSetting);
	  thermoLastUpdate = eventCount;
	  thermoCount = 0;
	}
    }

 
  updatePosition(particle); //update particle's position
  double factor = sqrt(temperature);
  //assign values for each component of the particle's velocity from Gaussian
  for(size_t dim (0); dim < 3; ++dim) 
      particle.v[dim] = RNG.var_normal() * factor;
  //cerr << "old v = (" << oldv.x << ", " << oldv.y << ", " << oldv.z << " ) "
  //     << "new v = (" << particle.v.x << ", " << particle.v.y << ", " << particle.v.z << ")" << endl;
  //cerr << "Change in kinetic energy of system " << 0.5 * particle.mass * (particle.v.dotProd(particle.v)-oldv.dotProd(oldv)) << endl;
  int particleNo = -1;
  double t_min_thermo = calcThermoTime(RNG, particleNo);
  if(particleNo == -1)
    {
      cerr << "calcThermoTime has not returned a valid particle No" << endl;
      exit(1);
    }
  events[numberParticles] = eventTimes(t_min_thermo, particleNo, -1, -1, eventTimes::THERMOSTAT);

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
  for(it_particle particle = particles.begin(); particle < particles.end(); ++particle) //loop through particles and calculate sum of x and y velocites
    for(size_t dim (0); dim < 3; ++dim)
      sum[dim] += particle->v[dim];
  
  for(it_particle particle = particles.begin(); particle < particles.end(); ++particle) //loop through particles and calculate sum of x and y velocites
    for(size_t dim (0); dim < 3; ++dim)
      particle->v[dim] -= sum[dim] / particles.size(); //reduce velocity by 1/number of particles of the total sum

}
