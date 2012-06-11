//----Includes----
#include <iostream>     //allows input/output to screen
#include <fstream>      //allows file io
#include <vector>       //allows use of vector structures
#include <string>       //allows use of strings
#include <time.h>       //allows use of timing functions
#include <unistd.h> // Header file for sleeping.
#include <stdlib.h> //allows exiting code
#include <math.h>
#include <GL/freeglut.h>    //allows use of the GLUT Library
#include <GL/gl.h>	//allows use of 3d gfx
#include <GL/glu.h>	//allows use of 3d gfx utilities

//----Structures----
struct Particle
{
  Particle(double xin, double yin, double zin): x(xin), y(yin), z(zin) {}
  double x;
  double y;
  double z;
};

struct Event
{
  Event(double time, std::vector<Particle>& p): t(time),particles(p) {}
  double t;
  std::vector<Particle> particles;
};
std::vector<Event> events;
#include "fileBuffer.h"
//----Defines----
#define ESCAPE 27 //ascii code for the escape code

//----Function Declarations----
void keyPress();
void drawCube(float, float, float, float);
