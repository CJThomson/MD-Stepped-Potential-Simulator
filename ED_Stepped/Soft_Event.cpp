//----Program Includes----
#include "Declares.h" //all includes and function declares for event sim
//Physical Properties:
double density = 0.85;
double temperature = 1.34; //temperature of the system
//Simulation:
int numberParticles = 256; //number of particles
const int numberEvents = 1.5e+6;
int eventCount = 0;
double length = pow(numberParticles/density, 1.0 / 3.0);
int number_of_runs = 1;
//const double length = 6.0;
CVector3 systemSize(length,length,length); //size of the system
double t = 0;
const bool initFile = false; //use an init file instead of random generated values
const bool overwriteInit = false; //create a new init file
std::vector<Steps> steps; //create a vector to store step propeties
const int noCells = 3;
int no_of_steps = 10;
double r_cutoff = 3.0;

//Thermostat:
bool thermostat = true; //use a thermostat
const size_t thermoFreq = 100; //frequency thermostat rate is updated
size_t thermoCount = 0; //counter to change thermostat rate
size_t thermoLastUpdate = 0;
double thermoMeanFreeTime = 0.0005;
//const int thermoOff =3.5e+6;
const int thermoOff = numberEvents;
double thermoSetting = 0.05;

//Reduced Unit Definitions
const double mass = 1; //mass of a particle
const double radius = 0.5; //radius of a particle (set for diameter = 1)
const double lj_sigma = 1.0;
const double lj_epsilon = 1.0;
double Stepper::lj_eps = 1.0;
double Stepper::lj_sig = 1.0;
double Stepper::beta = 1.0 / temperature;
Stepper::StepHeight height_type = Stepper::VIRIAL;
Stepper::StepWidth width_type = Stepper::EXPECTEDFORCE;
//Logging:
const int psteps = 50; //frequency of output to file
const int writeOutLog = 0;//level of outLog, 0 = nothing, 1 = event discriptions, 2 = full

//Measuring Properties
const int startSampling = 5e+5; //step number to start taking samples
const int sample_interval = 1;
const double rdf_interval = 0.05;
const int diff_interval = 20;
int readingsTaken = 0;
int rdfReadings = 0;
double startSampleTime = 0;
double currentK = 0;
double currentU = 0;
const int noBins = 600; //number of radial bins
//const double maxR  = 0.5 * std::min(systemSize.x, std::min(systemSize.y, systemSize.z)); //maximum radial distribution considered;
const double maxR = 3.0;
double rdf_d[noBins]; //radial distribution values
std::vector<Diffusion> coDiff; //coefficient of diffusion over

//----Time Averages----
double TA_rdf_d[noBins];
std::vector<std::pair<double, double> > icf;
std::vector<std::pair<double, double> > rdf_c;
double TA_v = 0;
double TA_U = 0;
double TA_T = 0;
double TA_p = 0;
double TA_v2 = 0;
double TA_U2 = 0;
double TA_T2 = 0;
double TA_p2 = 0;
double TA_tavg = 0;
double TA_tavg2 = 0;
Logger logger; //create an instance of the logger class

std::map<std::pair<int, int>, int> collStep;
using namespace std;
int main()
{
  string input;
  cout << " - - - Stepped Potential ED Simulator - - - " << endl;
  while(true)
    {
      cout << "input: ";
      cin >> input;
      if(input == "exit")
	exit(0);
      else if(input == "start")
	{
	  Stepper::beta = 1.0 / temperature;
	  Stepper stepper;
	  time_t startTime;
	  time(&startTime);
	  cout << "Starting simulation at " << asctime(localtime(&startTime)) << endl;
	  cout << "Generating " << no_of_steps << " Steps...";
	  //initSteps(); //step up system steps
	  stepper.generateSteps(no_of_steps, r_cutoff, height_type, width_type, steps);
	  logger.write_Steps(steps, temperature, density, numberParticles, height_type);
	  cout << " Complete" << endl;
	  vector<Results> results;
	  cout << "Running simulation " << number_of_runs << " times " << endl << endl;
	  for(size_t runs (0); runs < number_of_runs; ++runs)
	    {
	      cout << "\rRun number: " << runs + 1 << " of " << number_of_runs << endl;
	      
	      resetSim();
	      runSimulation(results, runs);
	    }
	  time_t endTime;
	  time(&endTime);
	  cout << "\rSimulation finished at " << asctime(localtime(&endTime)) 
	       << " and took " << round(difftime(endTime, startTime) / 60.0) << " mins" <<endl;
	  Results avgResults;
	  for(vector<Results>::iterator result = results.begin(); result != results.end(); ++result)
	    avgResults += *result;
	  avgResults *= 1.0 / number_of_runs;
	  cout << "\rAverage Results" << endl;
	  cout << "Temperature: " << avgResults.temperature << endl;
	  cout << "Pressure (discont): " << avgResults.pressure_d << endl;
	  cout << "Pressure (cont): " << avgResults.pressure_c << endl;
	  cout << "Potential Energy (discont): " << avgResults.potential_d << endl;
	  cout << "Potential Energy (cont): " << avgResults.potential_c << endl;
	  logger.write_Results(results, density, temperature, numberParticles);
	}
      else if(input == "temperature")
	cin >> temperature;
      else if(input == "density")
	cin >> density;
      else if(input == "particles")
	cin >> numberParticles;
      else if(input == "runs")
	cin >> number_of_runs;
      else if(input == "cutoff")
	cin >> r_cutoff;
      else if(input == "nosteps")
	cin >> no_of_steps;
      else
	cout << "Invalid input" << endl;
    }
}
void resetSim()
{
  eventCount = 0;
  length = pow(numberParticles/density, 1.0 / 3.0);
  systemSize = CVector3(length,length,length); //size of the system
  t = 0;
  thermostat = true; //use a thermostat
  thermoCount = 0; //counter to change thermostat rate
  thermoLastUpdate = 0;
  thermoMeanFreeTime = 0.0005;
  thermoSetting = 0.05;

  readingsTaken = 0;
  rdfReadings = 0;
  startSampleTime = 0;
  currentK = 0;
  currentU = 0;
  rdf_c.clear();
  icf.clear();
  for(size_t i(0); i < noBins; ++i)
    TA_rdf_d[i] = 0;
  TA_v = 0;
  TA_U = 0;
  TA_T = 0;
  TA_p = 0;
  TA_v2 = 0;
  TA_U2 = 0;
  TA_T2 = 0;
  TA_p2 = 0;
  TA_tavg = 0;
  TA_tavg2 = 0;
}

void runSimulation(vector<Results>& results, size_t runNumber)
{
  //variable declarations
  vector<CParticle> particles; //create a vector to store particle info
  vector<vector<eventTimes> > particleEL;
  vector<eventTimes> masterEL; //create a vector to store collision times
  vector<set<int> > neighbourList; //create a vector of sets to hold particles in  neighbour cell
  vector<set<int> > neighbourCell; //create a vector of sets to hold the cells neighbouring each cell
  cout << "Initialising Random Number Generators";
  CRandom RNG;
  RNG.seed();
  cout << "\rInitialising Particles ...";
  //Initialise the simulation
  if(initFile)
    initFromFile(particles);
  else
    initialise(particles, RNG ); // initialise the system
  numberParticles = particles.size();
  currentK = calcKinetic(particles);

  cout << "\rInitialising neighbour lists";
  int cells3 = pow(noCells, 3);
  neighbourList.resize(cells3);
  neighbourCell.resize(cells3);
  generateNeighbourCells(neighbourCell);
  generateNeighbourList(neighbourList, particles);

  cout << "\rInitialising Steps...";

  if(length * 0.5 < steps[steps.size()-1].step_radius)
    cout << "Warning system size less than cut-off radius" << endl;
  
  cout << "\rPopulating initial capture map";
  collStep.clear();
  for(it_particle p1 = particles.begin(); p1 != particles.end(); ++p1)
    for(it_particle p2 = p1 + 1; p2 != particles.end(); ++p2)
      calcStep(*p1,*p2);
  checkCaptureMap(particles);
  cout << "\rInitialising Output Logs...";
  //logger.initialise(Logger::LOCATIONS); //initialise location logger
  logger.initialise(Logger::OUTPUTLOG);


  cout << "\rGenerating event list...";
  particleEL.resize(particles.size()); 
  masterEL.resize(particles.size() + 2);//one event for each particle plus thermostat
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
  masterEL[particles.size()] = eventTimes(rdf_interval, -1, -1, -1, eventTimes::RDF);
  masterEL[particles.size() + 1] = eventTimes(t_min_thermo, particleNo, -1, -1, eventTimes::THERMOSTAT);
  
  logger.write_Location(particles, 0, systemSize); //write initial values to the log

  cout << "Starting Simulation with " << numberParticles << " particles with a density of "
       << density << " at a target temperature of " << temperature << " in a box with side length of " << length << endl;
  for(;eventCount < numberEvents;)
    {
      bool collEvent = false;
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
		if(writeOutLog >= 1)
		  logger.outLog << "Invalid Collision between particles: " << next_event.particle1 << " & " << next_event.particle2 << endl;
	      }
	    else //if a valid event
	      {
		collEvent = true;
		freeStream(dt);
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
		
		double temp_pFlux = calcVelocity(particles, next_event);		
		TA_v += temp_pFlux;
		TA_v2 += temp_pFlux * temp_pFlux;

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

	  freeStream(dt);
	  getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
	  break;
	case eventTimes::NEIGHBOURCELL:
	  {
	    freeStream(dt);
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
	    freeStream(dt);
	    if(writeOutLog >= 1)
	      logger.outLog << "Thermostat event for particle: " << next_event.particle1 << " at time = "<< t<< endl;
	    runThermostat(particles[next_event.particle1], RNG, masterEL);
	    ++particles[next_event.particle1].collNo;
	    getEvent(particles[next_event.particle1], particles, particleEL, masterEL, neighbourList, neighbourCell);
	    break;
	  }
	case eventTimes::RDF:
	  {
	    freeStream(dt);
	    if(eventCount > startSampling)
	      {
		if(writeOutLog >=1)
		  logger.outLog << "Measuring RDF at time " << t << endl;
		++rdfReadings;
		calcRadDist(particles);
		for(size_t i = 0; i < noBins; ++i)
		  TA_rdf_d[i] += rdf_d[i];
	      }
	    masterEL[particles.size()] = eventTimes(t + rdf_interval, -1, -1, -1, eventTimes::RDF);
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
	      TA_v = 0; 
	      TA_T = 0;
	      TA_U = 0;
	    }

	  if(eventCount % diff_interval == 0)
	    coDiff.push_back(Diffusion(calcDiff(particles, t),t));

	  if(eventCount % (int) ceil(numberEvents / 1000) == 0)
	    //if(eventCount % 1 == 0)
	    {
	      double E_temp = TA_T / (t - startSampleTime);
	      double E_pot = TA_U / (t - startSampleTime);
	      double E_press_ideal = density * E_temp;
	      double E_press_coll = mass * density * TA_v / (numberParticles * 3.0 * (t - startSampleTime));


	      cout << "\rN:" << eventCount << "/" << numberEvents 
		   << " T: " << currentK / (1.5 * numberParticles)
		   << " <T>: " << E_temp
		   << " U: " << currentU
		   << " <U>: " << E_pot
		   << " <P>: " << E_press_ideal + E_press_coll
		//<< " <P_coll>: " << E_press_coll
		// << " TE: " << currentK + currentU
		   << flush;
	      
	    }
	}
    }

  cout << "\rSimulation Complete";
  checkCaptureMap(particles);
  double deltaR = maxR / noBins; //width of each shell
  
  for(int i = 0; i < noBins; ++i)
    {
      double volShell = 4.0 / 3.0 * M_PI * (pow(deltaR * (i + 1), 3) - pow(deltaR * i, 3));
      TA_rdf_d[i] /= (0.5 * numberParticles * rdfReadings * volShell * density);
    }

  cout << "\rGenerating Continuous g(r), U, P";
  double E_temp = TA_T / (t - startSampleTime);
  indirectCorr(E_temp);
  continuousRDF(E_temp);
  double cont_P = continuousP(E_temp);
  double cont_U = continuousU();
  double E_pot = TA_U / (t - startSampleTime);
  double E_press = density * E_temp 
    + mass * density * TA_v / (numberParticles * 3.0 * (t - startSampleTime));
  results.push_back(Results(E_temp, E_press, cont_P, E_pot, cont_U));
  //Output time averages
  if(false)
    {
      cout << "Time Averages:" << endl;
      {
	// = Mean free time
	double freeTime = t * numberParticles / (2 * eventCount);
	double freeTime2 = TA_tavg2 * numberParticles / (2 * readingsTaken);
	double sd_freeTime = sqrt(fabs(freeTime2 - freeTime * freeTime));
	cout << "Mean free time: " << freeTime << " (" << sd_freeTime << ")" << endl;
	
	// = Potential Energy
	double E_pot = TA_U / (t - startSampleTime);
	double E_pot2 = TA_U2 / (t - startSampleTime);
	double sd_pot = sqrt(fabs(E_pot2 - E_pot * E_pot));
	cout << "Ud: " << E_pot << "(" << sd_pot << ")" << endl;
	cout << "Uc: " << cont_U << endl;
	
	// = Temperature
	double E_temp = TA_T / (t - startSampleTime);
	double E_temp2 = TA_T2 / (t - startSampleTime);
	double sd_temp = sqrt(fabs(E_temp2 - E_temp * E_temp));
	cout << "T: " << E_temp << "(" << sd_temp << ")" << endl;
	
	// = Momentum flux
	double E_mf = TA_v / readingsTaken;
	double E_mf2 = TA_v2 / readingsTaken;
	double sd_mf = sqrt(fabs(E_mf2 - E_mf * E_mf));
	cout << "<r.v> " << E_mf << "(" << sd_mf << ")" << endl;
	
	// = Pressure
	double E_press = density * E_temp 
	  + mass * density * TA_v / (numberParticles * 3.0 * (t - startSampleTime));
	double sd_press = sqrt( pow(density * sd_temp, 2) 
				+ pow(mass * density / (numberParticles * 3.0 * (t - startSampleTime))  * sd_mf, 2));
	
	cout << "Pd: " << E_press  << "(" << sd_press << ")" << endl;
	cout << "Pc: " << cont_P << endl;
      }
    }
  cout << "\rWriting Results...";
  logger.write_contRDF(rdf_c, density, temperature, numberParticles);
  logger.write_ICF(icf, density, temperature, numberParticles);
  logger.write_RadDist(TA_rdf_d, noBins, deltaR, density, temperature, numberParticles); //write radial distribution file
  logger.write_Diff(coDiff); //write coefficient of diffusion file
  logger.write_Location(particles, t, systemSize); //write final values to the log
  if(overwriteInit)
    logger.write_Init(particles); //write final values to log to allow
  cout << "\rAll Tasks Complete";
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

  return kinetic;
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
    rdf_d[i] = 0; //rezero array

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
	      ++rdf_d[index];
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
	  currentU += steps[i].step_energy;
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
	    CVector3 deltav1 = (A / mass) * r12.normalise();
	    particle[p1].v += deltav1;
	    particle[p2].v -= deltav1;

	    if(it_map == collStep.end()) //if no collision state found then particles must be outside outer step
	      collStep.insert(pair<pair<int, int>, int>(pair<int, int>(p1, p2), steps.size() - 1)); //insert pair into collStep map
	    else
	      --(it_map->second); //move particles in one step
	    currentU -= dU;
	    currentK += dU;
	    return r12.dotProd(deltav1);
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >= 2)
	      logger.outLog << "Well Bounce" << endl;
	    CVector3 deltav1 = - vdotr * r12.normalise();
	    particle[p1].v += deltav1;
	    particle[p2].v -= deltav1;
	    return r12.dotProd(deltav1);
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
	    CVector3 deltav1 = A / mass * r12.normalise();
	    particle[p1].v += deltav1;
	    particle[p2].v -= deltav1;
	    if(it_map->second == steps.size() -1) //if particle is leaving outermost step
	      collStep.erase(it_map);
	    else
	      ++(it_map->second); //move particles in one step
	    currentU -= dU;
	    currentK += dU;
	    return r12.dotProd(deltav1);
	  }
	else //if bounce occurs
	  {
	    if(writeOutLog >=2)
	      logger.outLog << "Well Bounce" << endl;
	    CVector3 deltav1 = - vdotr * r12.normalise();
	    particle[p1].v += deltav1;
	    particle[p2].v -= deltav1; 
	    return r12.dotProd(deltav1);
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
  steps.clear();
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

  /*steps.push_back(Steps(0.966872, 4.7605 )); //step 0
  steps.push_back(Steps(1.20474, -0.75868)); //step 1
  steps.push_back(Steps(1.38391, -0.670859)); //step 3
  steps.push_back(Steps(1.55901, -0.358529)); //step 4
  steps.push_back(Steps(1.7183, -0.19731)); //step 5
  steps.push_back(Steps(1.86003, -0.118788)); //step 6
  steps.push_back(Steps(1.98651, -0.0777275)); //step 7
  steps.push_back(Steps(2.10054, -0.0542997)); //step 8
  steps.push_back(Steps(2.20443, -0.0398837)); //step 9
  steps.push_back(Steps(2.3, -0.0304544)); //step 10*/

  /*steps.push_back(Steps(0.966872, 2.08615)); //step 0
  steps.push_back(Steps(1.20474, -0.813607)); //step 1
  steps.push_back(Steps(1.38391,   -0.675784)); //step 3
  steps.push_back(Steps(1.55901, -0.360147)); //step 4
  steps.push_back(Steps(1.7183, -0.197681)); //step 5
  steps.push_back(Steps(1.86003, -0.118882)); //step 6
  steps.push_back(Steps(1.98651, -0.0777557)); //step 7
  steps.push_back(Steps(2.10054,-0.0543097)); //step 8
  steps.push_back(Steps(2.20443, -0.0304544)); //step 9
  steps.push_back(Steps(2.3, -0.0304544)); //step 10*/


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
  currentK += 0.5 * particle.mass * (particle.v.dotProd()-oldv.dotProd());
  int particleNo = -1;
  double t_min_thermo = calcThermoTime(RNG, particleNo);
  if(particleNo == -1)
    {
      cerr << "calcThermoTime has not returned a valid particle No" << endl;
      exit(1);
    }
  events[numberParticles + 1] = eventTimes(t_min_thermo, particleNo, -1, -1, eventTimes::THERMOSTAT);

}
void freeStream(double dt)
{
  t += dt; //update system time
  double temp_temperature = currentK / (1.5 * numberParticles);
  TA_T += temp_temperature * dt;
  TA_U += currentU / numberParticles * dt;
  if(eventCount > startSampling)
    {
      if(startSampleTime <0) {cerr<< "ERROR" << endl; exit(1);}
      // = Mean Free Time
      TA_tavg2 += dt * dt;
      // = Temperature
      TA_T2 += temp_temperature * temp_temperature * dt;

      // = Potential
      TA_U2 += pow((currentU / numberParticles),2) * dt;
      ++readingsTaken;

    }
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
void checkCaptureMap(vector<CParticle> &particles)
{
  //CAPTURE TEST
  cout << "\rChecking the capture map ...";
  for(int i = 0; i < particles.size();++i)
    for(int j = i + 1; j < particles.size(); ++j)
      {
	CVector3 r12 = particles[i].r - particles[j].r;
	CVector3 v12 = particles[i].v - particles[j].v;
	applyBC(r12);
	map<pair<int, int>, int>::const_iterator it_map = collStep.find(pair<int,int>(i, j)); //find collision state of particles
	double distance = r12.length();
	int step = steps.size();
	for(int k = 0; k < steps.size(); ++k)
	  {
	    if(distance <= steps[k].step_radius)
	      { step = k; break; }
	  }
	
	{
	  if(it_map == collStep.end())
	    {
	      if (step < steps.size())
		cerr << "particles " << i << " and " << j << ", step=" << step << ", are not in capture map" << endl;
	    }
	  else
	    {
	      if (step == steps.size())
		cerr << "particles " << i << " and " << j << ",map_step=" << it_map->second << ", should not be in capture map" << endl;
	      else if (it_map->second != step)
		cerr << "particles " << i << " and " << j << ", step=" << step << ",map_step=" << it_map->second << ", are in wrong step" << endl;
	    }
	  
	  break;
	}

      }
}

void indirectCorr(double T)
{
  for(size_t i(0); i < noBins; ++i)
    {
      double distance = i * maxR / noBins;
      double potential_d = 0;
      double stepEnergy = 0;
      double stepRad = 0;
      for(size_t j(0); j < steps.size(); ++j)
	{
	  if(distance < steps[j].step_radius)
	    {
	      stepRad = steps[j].step_radius;
	      stepEnergy = steps[j].step_energy;
	      break;
	    }
	}
      if((distance < stepRad && (i+1) * maxR / noBins < stepRad) || stepRad == 0)
	{
	  //	  cerr << distance << " - " << stepRad << " - " << (i+1) * maxR / noBins << endl;
	  icf.push_back(pair<double, double> (distance, TA_rdf_d[i] * exp(stepEnergy / T)));
	}
    }
  icf[0] = pair<double, double>(0,0);

}
void continuousRDF(double T)
{
  for(vector<pair<double, double> >::iterator i = icf.begin(); i != icf.end(); ++i)
    {
      double distance = i->first;
      double potential_c = 4.0 * lj_epsilon * (pow(lj_sigma/distance, 12) - pow(lj_sigma/distance, 6)); 
      rdf_c.push_back(pair<double, double> (distance, i->second / exp(potential_c / T)));
    }
  rdf_c[0].second = 0;
}

double continuousU()
{
  //integration using trapezoidal rule
  double sum = 0;
  for(vector<pair<double, double> >::iterator i = rdf_c.begin(); i != rdf_c.end() - 1; ++i)
    {
      vector<pair<double, double> >::iterator j = i + 1;
      if(j == rdf_c.end())
	cerr<< "ERROR in continuous U" << endl;
      double u1 = 4.0 * lj_epsilon * (pow(lj_sigma/i->first, 12) - pow(lj_sigma/i->first, 6)); 
      double u2 = 4.0 * lj_epsilon * (pow(lj_sigma/j->first, 12) - pow(lj_sigma/j->first, 6)); 
      if(i->first == 0)
	u1 = 0;
      double f1 = u1 * i->second * i->first * i->first;
      double f2 = u2 * j->second * j->first * j->first;
      sum += (j->first - i->first) * (f1 + f2);
    }
  return M_PI * density * sum;
}

double continuousP(double T)
{
  double sum = 0;
  for(vector<pair<double, double> >::iterator i = rdf_c.begin(); i != rdf_c.end() - 1; ++i)
    {
      vector<pair<double, double> >::iterator j = i + 1;
      if(j == rdf_c.end())
	cerr<< "ERROR in continuous P" << endl;
      double u1 = 24.0 * lj_epsilon / lj_sigma * (2.0 * pow(lj_sigma/i->first, 13) - pow(lj_sigma/i->first, 7)); 
      if(i->first == 0)
	u1 = 0;
      double u2 = 24.0 * lj_epsilon / lj_sigma * (2.0 * pow(lj_sigma/j->first, 13) - pow(lj_sigma/j->first, 7)); 

      double f1 = u1 * i->second * i->first * i->first * i->first;
      double f2 = u2 * j->second * j->first * j->first * i->first;
      sum += (j->first - i->first) * (f1 + f2);
    }
  return density * T + 1.0 / 3.0 * M_PI * density * density * sum;
}
