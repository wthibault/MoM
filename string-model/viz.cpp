//
// viz.cpp - realtime vibrating string experimental apparatus
// 
/*
 * trackball.cpp
 * Dave Rogers, modified from
 * Tebo: trackball.py
 * 14-Dec-2010 - modified by Bill Thibault to draw a FDP'd graph
 */


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <vector>
#include <omp.h>
#include <GL/glut.h>

#include "RtAudio.h"

#include "util.h"

using namespace std;

//**** GLOBALS ******//

// WINDOW SIZE
GLint g_win_width = 800;
GLint g_win_height = 600;

GLfloat g_start_point[3];
GLfloat g_trackball_transform[16];

string g_usage;

bool  autorotate = false;
float autorotateAngle = 0.0;

class StringModel;
StringModel               *theString;
RtAudio::StreamParameters parameters;
unsigned int              sampleRate = 44100;
unsigned int              bufferFrames = 256; // 256 sample frames
RtAudio dac;
//
//////

//
// physical string model - 
//


struct StringModel {
  StringModel ( int n, float _Ktension, float _Kdamping, int _stepspersample )
    : numMasses(n), 
      Ktension(_Ktension),
      Kdamping(_Kdamping),
      simulationStepsPerSample(_stepspersample),
      vibratorOn ( false ),
      vibratorFreq ( 100.0f ),
      vibratorAmplitude ( 0.001f ),
      vibratorPhase ( 0.0 )
      
  {
    std::cout << "StringModel N,K_t,K_d,ss: " 
	      << numMasses << ',' 
	      << Ktension  << ',' 
	      << Kdamping  << ',' 
	      << simulationStepsPerSample
	      << std::endl;

    // allocate displacement and velocity arrays

    y = new float[numMasses];
    yold = new float[numMasses];
    v = new float[numMasses];

    // initialize displacements and velocities

    for (int i = 0; i < numMasses; i++ ) {
      v[i]  = 0.0f;
      yold[i] = y[i] = 0.0f;
    }

    // init the mutex
    pthread_mutex_init ( &lock, NULL );

    // init the RNG
    seed = (unsigned int) time(NULL);

    // pluck it
    pluck();
  }

  ~StringModel() { delete y; delete yold, delete v; }


  void print() {
    pthread_mutex_lock ( &lock );
    std::cout << "N,K_t,K_d,ss: " 
	      << numMasses << Ktension << Kdamping << simulationStepsPerSample
	      << std::endl;
    for (int i = 0; i < numMasses; i++) {
      std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    pthread_mutex_unlock ( &lock );
  }


  void reset() {
    pthread_mutex_lock ( &lock );
    for (int i = 0 ; i < numMasses; i++ ) {
      v[i] = y[i] = yold[i] = 0.0;
    }
    pthread_mutex_unlock ( &lock );
  }


  void pluck() {
    int pluckAt = float(rand_r(&seed) / float(RAND_MAX)) * (numMasses-2) + 1;
    //    std::cout << pluckAt << std::endl;
    pthread_mutex_lock( &lock );
    float maxDisp = 0.001;
    float upSlope = maxDisp / pluckAt;
    float downSlope = maxDisp / (numMasses-pluckAt);
    for (int i = 1; i< numMasses-2;i++ ) {
      if (i <= pluckAt )
	yold[i] = i*upSlope;
      else
	yold[i] = maxDisp - (i-pluckAt)*downSlope; 
    }
    pthread_mutex_unlock ( &lock );
  }

  void pluckvel() {
    int pluckAt = float(rand_r(&seed) / float(RAND_MAX)) * (numMasses-2) + 1;
    //    std::cout << pluckAt << std::endl;
    pthread_mutex_lock( &lock );
    v[pluckAt] = 3e-4 * Ktension;
    pthread_mutex_unlock ( &lock );
  }

  void toggleVibrator() {
    vibratorOn = !vibratorOn;
  }

  static int audioCallback ( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
			   double streamTime, RtAudioStreamStatus status, void *userData );

  void computeSamples ( float *outputBuffer, unsigned int nBufferFrames );
  int simulationStepsPerSample;
  int numMasses;
  float Ktension;
  float Kdamping;
  float *y, *yold, *v;
  unsigned int seed;
  pthread_mutex_t lock;
  bool   vibratorOn;
  float vibratorFreq;
  float vibratorAmplitude;
  float vibratorPhase;
  float t;
};

// pure c++ audio computation
class Parameter {
  Parameter() : targetValue(0.0f),
		previousValue ( 0.0f )
  {}
  inline float getValueAt(float u) {
    return lastValue = u * targetValue + (1.0f-u) * previousValue;
  }
  inline void setTargetValue( float target ) {
    targetValue = target;
    previousValue = lastValue;
  }
  inline void setCurrentValue ( float v ) {
    previousValue = lastValue = v;
  }
  float targetValue;
  float previousValue;
  float lastValue;
};

void
StringModel::computeSamples ( float *soundout, unsigned int nBufferFrames )
{
  int iters = nBufferFrames * simulationStepsPerSample;

  for ( int t = 0; t < iters; t++ ) {

    float sum = 0;
    int n = numMasses;
    int i;
    float accel;
    if ( vibratorOn ) {
      // XXX interpolate params
      y[0] = vibratorAmplitude * sin ( vibratorFreq * vibratorPhase );
      vibratorPhase += (1.0 / sampleRate) / simulationStepsPerSample;
      while ( vibratorPhase > 2*M_PI )
	vibratorPhase -= 2*M_PI;
    }

    // XXX not working : #pragma omp parallel for reduction(+:sum) private(accel)
    for ( i = 1; i < n-1; i++ ) {
	accel = Ktension * (yold[i+1] + yold[i-1] - 2*yold[i]);
	v[i] += accel;
	v[i] *= Kdamping;
	y[i] = yold[i] + v[i];
	sum = sum + y[i];
    }

    float *tmp = y;
    y = yold;
    yold = tmp;

    if ( t % simulationStepsPerSample == 0 ) {
      if (fabs(sum) > 1.0) {
	sum = sum < 0.0 ? -1.0 : 1.0;
	//	std::cout << "! " << sum << std::endl;
      }
#if 1
      // summed output
      *soundout++ = sum;
      *soundout++ = sum;
#else
      // per channel pickup placement
      *soundout++ = 10.0f * yold[numMasses/13];
      *soundout++ = 10.0f * yold[numMasses/5];
#endif
    }

  }

}

// the RtAudio callback

int
StringModel::audioCallback ( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	 double streamTime, RtAudioStreamStatus status, void *userData )
{

  if (status == RTAUDIO_INPUT_OVERFLOW)
    std::cout << "string -- overflow" << std::endl;

  float *soundout = (float *)outputBuffer;
  StringModel *s = (StringModel *)userData;

  pthread_mutex_lock(&(s->lock));

  s->computeSamples ( static_cast<float*>(outputBuffer), nBufferFrames );
  
  pthread_mutex_unlock ( &(s->lock) );
  return 0;
}

//////////////////////////////
// 


// /////////////////////////////
// // each worker thread gets a portion of the string
// pthread_cond_t  workReady = PTHREAD_COND_INITIALIZER;
// pthread_cond_t  workDone = PTHREAD_COND_INITIALIZER;
// float *audioBuffer;
// int nBufferFrames;
// int numThreads;
// void *threadFunc ( void * arg )
// {
//   long int id = (long int)(arg);
//   int chunkSize = theString->numMasses / numThreads;
//   int myStart = id * chunkSize;

//   if (myStart + chunkSize > theString->numMasses )
//     chunkSize = theString->numMasses - myStart;

//   // wait for data ready

//   /* for each buffer frame:
//        render my section into y
//        compute my sum  ( master does the global sum )
       
//   */
// }




//////


void init(int,char**);
void map_to_sphere(GLfloat out[3], int x, int y);
void mouse(int button, int state, int x, int y);
void mouse_motion(int x, int y);
void draw_axes();
void display();
void reshape(int, int);
void keyboard(unsigned char, int, int);
int main(int, char **);

void init(int argc, char **argv) {

  /////
  theString = new StringModel ( 1000, 0.5, 0.99999, 8 );

  // *** test the rtaudio callback
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
    
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  sampleRate = 44100;
  bufferFrames = 256; // 256 sample frames

  try { 
    dac.openStream ( &parameters, 
		     NULL, 
		     RTAUDIO_FLOAT32,
		     sampleRate, 
		     &bufferFrames, 
		     // &stringmodel, 
		     StringModel::audioCallback,
		     (void *)theString );
    dac.startStream();
  }
  catch ( RtError& e ) {
    std::cout << "\nexception on dac:\n";
    e.printMessage();
    exit(0);
  }

  //////
    
  GLfloat pos[] = {5.0, 5.0, 10.0, 0.0};
  glLightfv(GL_LIGHT0, GL_POSITION, pos);
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_FOG);
  float FogCol[3]={0.0,0.0,0.0};
  glFogfv(GL_FOG_COLOR,FogCol); 
  glFogi(GL_FOG_MODE, GL_LINEAR);
  glFogf(GL_FOG_START, 10.0f);
  glFogf(GL_FOG_END, 40.f);

  glClearColor (0.0, 0.0, 0.0, 0.0);

  set_to_ident(g_trackball_transform);

  std::cout << "\nPlaying ... press \n";
  std::cout << "t to tighten\n";
  std::cout << "l to loosen\n";
  std::cout << "p to pluck\n";
  std::cout << "r to reset\n";
  std::cout << "d to dump velocities\n";
  std::cout << "ESC to quit.\n";
}

/**
 * Map mouse coords x, y to Point on sphere
 */
void map_to_sphere(GLfloat out[3], int x, int y) {
	GLfloat _x, _y;
	_x = 2 * (x/(GLfloat)g_win_width) - 1;
	_y = 2 * ((g_win_height-y)/(GLfloat)g_win_height) - 1;

	GLfloat p[3] = {_x, _y, 0};
	GLfloat len = magnitude(p);
	GLfloat v[3];
	if(len > 1) {
		copypv(v, p);
		normalize(v);
	} else {
		v[0] = _x;
		v[1] = _y;
		v[2] = sqrt(1 -len*len);
//		v[2] = 1 -len*len;
	}
	copypv(out, v);
}

void mouse(int button, int state, int x, int y) {
	if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		map_to_sphere(g_start_point, x,y );
	}
	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
	  float amt = g_win_height - y;
	  amt *= 10.0;
	  theString->vibratorFreq = amt;
	}
}

void mouse_motion(int x, int y) {
	GLfloat new_point[3];
	map_to_sphere(new_point, x, y);

	GLfloat axis[3];
	cross(axis, g_start_point, new_point);
	copypv(g_start_point, new_point);

	GLfloat len = magnitude(axis);
	normalize(axis);
	GLfloat angle = ((GLfloat)asin(len));
	angle *= 180/PI;

	if(angle > .00001) {
		GLfloat transform[16];
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glRotatef(angle, axis[0], axis[1], axis[2]);
		glGetFloatv(GL_MODELVIEW_MATRIX, transform);

		glPopMatrix();

		GLfloat out[16];
		mult_matrixf(out, g_trackball_transform, transform);
		copym(g_trackball_transform, out);

		glutPostRedisplay();
	}
}



void draw_axes() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glDisable(GL_LIGHTING);
	glColor3f(1,0,0);
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(1,0,0);
	glEnd();
	glColor3f(0,1,0);
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(0,1,0);
	glEnd();
	glColor3f(0,0,1);
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(0,0,1);
	glEnd();
	glPopAttrib();

}


void drawString (StringModel *s)
{
  int n = s->numMasses;
  float *y = new float[n];
  float xscale = 1.0 / n;
  float xoffset = -0.5;
  float yscale = 20.0;
  float yoffset = 0;//-500;
  float velscale = 10000;

  pthread_mutex_lock(&s->lock);
  for (int i = 0; i < n; i++ )
    y[i] = s->y[i];
  pthread_mutex_unlock(&s->lock);

  glColor3f ( 1,1,1 );
  glLineWidth ( 3.0 );
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < n; i++ ) {
    float vel = fabs(s->v[i]);
    //    glColor3f ( velscale * vel, velscale*2 * vel, 1.0 );
    glVertex3f ( i*xscale + xoffset, 0, y[i]*yscale + yoffset );
  }  
  glEnd();

  if ( s->vibratorOn ) { 
    glColor3f ( 1, 0, 0 );
    glLineWidth ( 7.0 );
    glBegin ( GL_LINES );
    glVertex3f ( -10*xscale + xoffset, 0, y[0]*yscale + yoffset );
    glVertex3f ( xoffset, 0, y[0]*yscale + yoffset );
    glEnd ();
  }

  delete y;

}


void display(void)
{
   //
   // XXX TODO - move to OpenGL 3.2 or above (glDrawArrays(...))
   //

   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

   glPushMatrix();
   glMultMatrixf(g_trackball_transform);

   //   draw_axes();

   if (autorotate) {
     autorotateAngle += 0.05;
   }
   glRotatef(autorotateAngle, 0,1,0 );

   glRotatef(90, 1, 0, 0);

   /////

   drawString ( theString );

   /////


   glPopMatrix();

   glutSwapBuffers();
}

void reshape (int w, int h) {
  g_win_width = w;
  g_win_height = h;
  glViewport (0, 0, (GLsizei) g_win_width, (GLsizei) g_win_height);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective(65.0, (GLfloat) w/(GLfloat) h, 0.1, 200.0);
  //glOrtho( 0, theString->numMasses, -0.01, 0.01, -10, 10);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef (0.0, 0.0, -1.0);
}

void keyboard (unsigned char key, int x, int y)
{
  switch (key) {
  case 'a':
    autorotate = ! autorotate;
    break;

  case 't' : theString->Ktension *= 1.05946; break;
  case 'l' : theString->Ktension *= 0.943876; break;
  case 'P' : theString->pluck(); break;
  case 'p' : theString->pluckvel(); break;
  case 'r' : theString->reset(); break;
  case 'd' : theString->print(); break;
  case 'v' : theString->toggleVibrator(); break;
  case 27: /* ESC */
    try {
      // Stop the stream
      dac.stopStream();
    }
    catch (RtError& e) {
      e.printMessage();
    }
	  
    if (dac.isStreamOpen()) dac.closeStream();
    exit(0);
    break;
  default: break;
  }
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitWindowSize (g_win_width, g_win_height);
	glutInitWindowPosition (320,240);
	glutCreateWindow ("String Model");
	init (argc,argv);
	cout << g_usage;

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc (mouse);
	glutMotionFunc(mouse_motion);

	glutIdleFunc(display);
	glutMainLoop();
	return 0;
}
