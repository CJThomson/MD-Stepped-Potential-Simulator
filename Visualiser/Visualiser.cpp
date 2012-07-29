#include "Declares.h"
using namespace std;

int window; //GLUT window

float camX = 0, camY = 0, camZ = 0; //camera position
float rotX = 0, rotY = 0; // rotate camera
int scrnHeight = 768; //height of window
int scrnWidth = 768; //width of window
float length = 0; //length of boundary box
int n = 0; //event number
int buffercount = 0;
int framerate = 400; //rate events are updated
bool bPause = true;
bool showall = false;
string title = "Particle Visualiser"; //window title
clock_t time_now, update_timer = 0;
GLUquadricObj *quadratic; //allows the drawing of spheres

fileBuffer fileBuf;
void InitGL(int Width, int Height) //initialise window
{
  GLfloat light_position0[] = {1.0f, 2.0f, 3.0f, 0.0f};
  GLfloat sphereColour[] = {0.5f, 0.5f, 1.0f, 1.0f};
  GLfloat sphereLight[] = {0.5f, 0.5f, 1.0f, 0.4f};
  glClearColor(0.7f, 0.7f, 1.0f, 0.0f); //set background colour
  glClearDepth(1.0); //allows clearing of the depth buffer
  glDepthFunc(GL_LESS); //depth testing buffer
  glLightfv(GL_LIGHT0, GL_POSITION,light_position0);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, sphereColour);
  glLightfv(GL_LIGHT1, GL_AMBIENT, sphereColour);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST); //enable depth testing
  glShadeModel(GL_SMOOTH); //enable smooth colour shading

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity(); //reset projection matrix
  //calc aspect ratio
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW); //reset modelview matrix
}

void ReSizeGLScene(int Width, int Height) //resize the window
{
  if (Height==0) //set height to 1 if height is 0, prevents divide by 0
    Height=1;

  glViewport(0, 0, Width, Height); //reset viewport, projection and modelview matrices

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,1000.0f);
  glMatrixMode(GL_MODELVIEW);
}

void DrawGLScene() //draw the scene
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //clear screen and depth buffer
  glLoadIdentity();
  glTranslatef(camX,camY,camZ); //move the camera
  //rotate camera
  glRotatef(rotX,1.0f,0.0f,0.0f);
  glRotatef(rotY,0.0f,1.0f,0.0f);
  glColor3f(0.0f,1.0f,0.0f); //set colour of the boundary box
  glPolygonMode( GL_FRONT_AND_BACK, GL_LINE ); //use wireframe models
  drawCube(-length * 0.5, -length * 0.5, -length * 0.5, length); //draw the boundary box
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  time_now = clock(); //take current time
  for(int x = -1; x < 2; ++x)
    for(int y = -1; y < 2; ++y)
      for(int z = -1; z < 2; ++z)
	{
	  if(showall)
	    {
	      glPushMatrix(); //push 2
	      glTranslated(x * length, y * length, z * length); //move sphere
	    }
	  for(int i = 0; i < events[buffercount].particles.size(); ++i)
	    {
	      glPushMatrix(); //push 3
	      glTranslated(events[buffercount].particles[i].x,events[buffercount].particles[i].y,events[buffercount].particles[i].z); //move sphere
	      gluSphere(quadratic,0.5f,20,20); //draw particle
	      glTranslated(-events[buffercount].particles[i].x,-events[buffercount].particles[i].y,-events[buffercount].particles[i].z); //move sphere
	      glPopMatrix(); //pop 3

	    }
	  if(showall)
	    glPopMatrix(); //pop 2
	}
  if(time_now - update_timer > CLOCKS_PER_SEC/framerate) //limit max framerate
    {
      //cout << time_now << " - " << update_timer << "c/s: " << CLOCKS_PER_SEC<< "frames: " << ceil((double)(time_now - update_timer)/CLOCKS_PER_SEC*framerate) <<endl;
      if(!bPause) //if simulation is not paused update scene
	{
	  int noFrames = ceil((double)(time_now - update_timer)/CLOCKS_PER_SEC*framerate);
    if(buffercount < 9)
      {
        fileBuf.bufferFile(10);
        buffercount = 0;
      }
    if(noFrames >= 10)
      noFrames = 1;
	  if(buffercount + noFrames < fileBuf.endFile) //if the last scene
	    buffercount += noFrames;
	  else if(buffercount < fileBuf.endFile)
	    ++buffercount;
	  else
	    glutLeaveMainLoop();
	}
      update_timer = time_now;
    }

  glutSwapBuffers(); //switch buffer (double buffering)
}
void keyUpNormal(unsigned char key, int x, int y)
{
  if(key == 'p')
    bPause = !bPause;
  if(key == 'j')
    showall = !showall;
}
void keyPressNormal(unsigned char key, int x, int y) //allows keypresses
{
  //usleep(100); //prevents hanging
  cout << key << endl;
  if(key == 27)
    {
      glutDestroyWindow(window);
      exit(0);
    }
  if(key == 'w')
    camY -= 0.2f;
  if(key == 's')
    camY += 0.2f;
  if(key == 'a')
    camX += 0.2f;
  if(key == 'd')
    camX -= 0.2f;
  if(key == 'q')
    camZ += 0.5f;
  if(key == 'e')
    camZ -= 0.5f;
  if(key == 'n' && buffercount<fileBuf.endFile-1)
    ++buffercount;
  if(key == 'b' && n>0)
    --n;
}
void keyPressSpecial(int key, int x, int y)
{
  switch(key)
    {
    case GLUT_KEY_UP:
      rotX += 1.5f;
      break;
    case GLUT_KEY_DOWN:
      rotX -= 1.5f;
      break;
    case GLUT_KEY_LEFT:
      rotY += 1.5f;
      break;
    case GLUT_KEY_RIGHT:
      rotY -= 1.5f;
      break;
    }
}

int main(int argc, char **argv)
{
  string input;
  cout << " - - :Particle Visualiser - - " << endl;
  while(true)
    {
      cout << "input: ";
      cin >> input;
      if(input == "exit")
        exit(0);
      else if(input == "framerate")
      {
        cout << "Enter new framerate: ";
        cin >> framerate;
      }
      else if(input == "showall")
      {
        cout << "Show all image (1/0): ";
        cin >> showall;
      }
      else if(input == "run")
      {
        n = 0; //restart the simulation from the beginning

        fileBuf.openFile();
        fileBuf.bufferFile(10);
        length = fileBuf.length;
        camZ = -5.0;
        quadratic = gluNewQuadric();          // Create A Pointer To The Quadric Object
        gluQuadricNormals(quadratic, GLU_SMOOTH);   // Create Smooth Normals
        glutInit(&argc, argv); //Initialise GLUT
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH); //set display mode

        glutInitWindowSize(scrnHeight, scrnWidth); //intialise window
        glutInitWindowPosition(50, 50); //set initial position

        window = glutCreateWindow(title.c_str()); //open window
        glutDisplayFunc(&DrawGLScene); //draw the scene
        glutIdleFunc(&DrawGLScene); //if idle just draw the scene

        glutReshapeFunc(&ReSizeGLScene); //if resized call our resize function
        glutIgnoreKeyRepeat(GLUT_KEY_REPEAT_OFF);
        glutSpecialFunc(&keyPressSpecial);
        glutKeyboardFunc(&keyPressNormal); //if keypress call our keypress function
        glutKeyboardUpFunc((&keyUpNormal));
        InitGL(scrnHeight, scrnWidth); //initailise window
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE,
          GLUT_ACTION_CONTINUE_EXECUTION); //on loop end don't close window
        glutMainLoop(); //main processing loop
      }
    }
  return 1;
}


void drawCube(float x, float y, float z, float sideLength)
{
  glBegin(GL_QUADS); //start drawing some quads
  //bottom
  glVertex3f(x, y, z);
  glVertex3f(x + sideLength, y, z);
  glVertex3f(x + sideLength, y, z + sideLength);
  glVertex3f(x, y, z + sideLength);
  //top
  glVertex3f(x, y + sideLength, z);
  glVertex3f(x + sideLength, y + sideLength, z);
  glVertex3f(x + sideLength, y + sideLength, z + sideLength);
  glVertex3f(x, y + sideLength, z + sideLength);
  //left
  glVertex3f(x, y, z);
  glVertex3f(x, y + sideLength, z);
  glVertex3f(x, y + sideLength, z + sideLength);
  glVertex3f(x, y, z + sideLength);
  //right
  glVertex3f(x + sideLength, y, z);
  glVertex3f(x + sideLength, y + sideLength, z);
  glVertex3f(x + sideLength, y + sideLength, z + sideLength);
  glVertex3f(x + sideLength, y, z + sideLength);
  //front
  glVertex3f(x, y, z);
  glVertex3f(x + sideLength, y, z);
  glVertex3f(x + sideLength, y + sideLength, z);
  glVertex3f(x, y + sideLength, z);
  //back
  glVertex3f(x, y, z + sideLength);
  glVertex3f(x + sideLength, y, z + sideLength);
  glVertex3f(x + sideLength, y + sideLength, z + sideLength);
  glVertex3f(x, y + sideLength, z + sideLength);
  glEnd();
}

