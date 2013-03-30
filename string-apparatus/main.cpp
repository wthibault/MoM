//
// main.cpp - main for string apparatus simulator with fltk-based gui
// // supercedes string_apparatus.cpp, kept as "pure glut" version is needed
//


#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef __APPLE__
#include <GL/glew.h>
#endif

#include "ssg.h"
#include "glm/glm.hpp"
#include "Camera.h"
#include "ParticleSystem.h"
#include "Trackball.h"
#include "StringModel.h"
#include "StringModelPrimitive.h"
#include "FFTPrimitive.h"
#include <fftw3.h>
#include "gui.h"

#ifdef USE_OSC
#include "osc.h"
#include <pthread.h>
#include <map>
#endif

using namespace glm;

//
// GLOBALS
//

// scene graph 
ModelNode *root;          

// physics model of string
StringModel *theString;

// audio output handle
RtAudio dac;

const char *shader = "PhongShading";

// window and widget constants
const int winWidth = 1024;
const int winHeight = 768;
const int offsetWidgets = winHeight / 2;
const int coarsefineHeight = (winHeight - offsetWidgets) / 7;

// audio params
const int sampleRate = 44100;
unsigned int bufferFrames = 256; // 256 sample frames ~ 5ms 
//unsigned int bufferFrames = 1024; // 1024 sample frames ~ 20ms 

// simulation params

const int simStepsPerSample = 2;
double initHangerMass = 20.0; // g
double initMassDensity = 0.1;  // g/m
double initDecayTime = 0.5;
double initVibFreq = 50.0;
double initVibAmp = 0.001;

#ifndef NUM_MASSES
#define NUM_MASSES 1000
#endif

///////////////////////////////////////////////////////////////////////



void 
keyboard (unsigned char key, int x, int y)
{
  float upCoarse = 1.06;
  float upFine = 1.001;
  float downCoarse = 1 / upCoarse;
  float downFine = 1 / upFine;
  float dampAdjustment = 1e-5;

  switch (key) {

  case 't' : 
    // coarse tension control - up
    theString->Ktension *= upCoarse;
    theString->Ktension = fmin ( 1.0, theString->Ktension );
    break;
  case 'T' : 
    // fine tension control - up
    theString->Ktension *= upFine;
    theString->Ktension = fmin ( 1.0, theString->Ktension );
    break;
  case 'g' : 
    // coarse tension control - down
    theString->Ktension *= downCoarse;
    break;
  case 'G' : 
    // fine tension control - down
    theString->Ktension *= downFine;
    break;


  case 'e' : 
    // coarse damping control - up
    //    theString->Kdamping *= upCoarse * dampAdjustment;
    break;
  case 'E' : 
    // fine damping control - up
    //    theString->Kdamping *= upFine / dampAdjustment;
    break;
  case 'd' : 
    // coarse damping control - down
    //    theString->Kdamping *= downCoarse / dampAdjustment;
    break;
  case 'D' : 
    // fine damping control - down
    //    theString->Kdamping *= downFine / dampAdjustment;
    break;



  case 'u' : 
    // coarse vibrator freq - up
    theString->vibratorFreq *= upCoarse;
    //    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;
  case 'U' : 
    // fine vibrator freq control - up
    theString->vibratorFreq *= upFine;
    //    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;
  case 'j' : 
    // coarse vibrator freq - down
    theString->vibratorFreq *= downCoarse;
    //    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;
  case 'J' : 
    // fine vibrator freq control - down
    theString->vibratorFreq *= downFine;
    //    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;


  case 'i' : 
    // coarse vibrator amp - up
    theString->vibratorAmplitude *= upCoarse;
    break;
  case 'I' : 
    // fine vibrator freq control - up
    theString->vibratorAmplitude *= upFine;
    break;
  case 'k' : 
    // coarse vibrator freq - down
    theString->vibratorAmplitude *= downCoarse;
    break;
  case 'K' : 
    // fine vibrator freq control - down
    theString->vibratorAmplitude *= downFine;
    break;



  case 'v' : theString->toggleVibrator(); break;
  case 'P' : theString->pluck(); break;
  case 'p' : theString->pluckvel(); break;
  case 'r' : theString->reset(); break;
    //  case 'd' : theString->print(); break;

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

  default:
    break;
  }
}


void
initOpacityValues ( float values[], int size )
{
  const float EPSILON = 1e-6;
  for ( int i = 0; i < size; i++ ) {
    float x = float(i) / float(size);
    values[i] = std::min ( 1.0, 1.0 / (M_PI * sqrt ( 1 - (2*x-1) * (2*x-1) )) );
  }
}


/////////////////////////////////////////////////////


class MyWindow : public Fl_Gl_Window {
  void draw();
  int handle(int);
  ModelNode *root;
  Camera     camera;
  Instance *loadBackgroundScene();
#ifdef USE_OSC
  OscParams          oscParams;
  pthread_mutex_t    oscMutex;
  pthread_t          oscThreadId;
  void               oscCollectInput();
#endif
public:
  MyWindow ( int x, int y , int w, int h, const char *L );
  void init();

};

MyWindow::MyWindow ( int x, int y , int w, int h, const char *L )
  : Fl_Gl_Window ( x, y, w, h, L ) , camera (), root ( NULL ) 
{
  mode ( FL_RGB8 | FL_DOUBLE | FL_ALPHA | FL_DEPTH | FL_MULTISAMPLE );
  camera.enableTrackball(true);
  camera.setDistance ( 2 );
  end();
#ifdef USE_OSC
  pthread_mutex_init ( &oscMutex, NULL );
#endif
}


Instance *
MyWindow::loadBackgroundScene()
{
  // be sure to get these from blender/string-apparatus.blend as the dimensions have to match numerous magic numbers in the gui
  std::vector<const char*> parts;
  parts.push_back( "objfiles/vibrator.obj" );
  parts.push_back( "objfiles/Desk.obj" );
  parts.push_back( "objfiles/MassHolder.obj" );
  parts.push_back( "objfiles/Cylinder.obj" );
  parts.push_back( "objfiles/VibratorStand.obj" );

  Instance *instance = new Instance();
  instance->setMatrix ( mat4() );

  Primitive *prim;
  //  prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );
  //prim = new ObjFilePrimitive ( "objfiles/string-scene.obj" );
  for (int i = 0; i < parts.size(); i++ ) {
    prim = new ObjFilePrimitive ( parts[i] );
    instance->addChild ( prim );
  }

  return instance;
}


void
MyWindow::init()
{
#ifndef __APPLE__
  glewExperimental=1;
  GLenum err = glewInit();
  if ( err != GLEW_OK ) {
    std::cout << "glewInit failed, aborting" << std::endl;
    exit(1);
  }
#endif


  Instance *instance = loadBackgroundScene();

  // the lights are global for all objects in the scene
  RenderingEnvironment::getInstance().lightPosition = vec4 ( 5,5,10,1 );
  RenderingEnvironment::getInstance().lightColor = vec4 ( 1,1,1,1 );
  // gravity center is in world coords, used by all ParticleSystem instances
  ParticleSystem::gravityCenter = vec3 ( 0,-1.9,0 );

  // create a material to use
  Material *mat = new Material;
  mat->ambient = vec4 ( 0.1, 0.1, 0.2, 1.0 );
  mat->diffuse = vec4 ( 0.5, 0.5, 0.1, 1.0 );
  mat->specular = vec4 ( 0.1, 0.1, 0.1, 1.0 );
  mat->shininess = 133.0;

  mat->program = mat->loadShaders ( shader );

  // set the instance as the scene root
  root = instance;

  // misc OpenGL state
  glClearColor (0.2, 0.2, 0.5, 1.0);
  glEnable(GL_DEPTH_TEST);
  glPointSize(1.0);
  glLineWidth(1.0);
  // for transparency
  glEnable ( GL_BLEND );
  glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

  //  glEnable ( GL_LINE_SMOOTH );
  //  glHint ( GL_LINE_SMOOTH, GL_NICEST );

  // the simulation 
  theString = new StringModel ( NUM_MASSES, 
				initHangerMass * 9.8,
				initMassDensity, 
				initDecayTime, 
				2.0,
				simStepsPerSample, 
				sampleRate, 
				bufferFrames );

  // the audio
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
    
  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  RtAudio::StreamOptions options;
  options.flags |= RTAUDIO_SCHEDULE_REALTIME;
  //  options.flags |= RTAUDIO_HOG_DEVICE;

  dac.showWarnings ( true );
  try { 
    dac.openStream ( &parameters, 
		     NULL, 
		     RTAUDIO_FLOAT64,
		     sampleRate, 
		     &bufferFrames, 
		     StringModel::audioCallback,
		     (void *)theString,
		     &options );
    dac.startStream();
  }
  catch ( RtError& e ) {
    std::cout << "\nexception on dac:\n";
    e.printMessage();
    exit(0);
  }

  // the rendering
  StringModelHistogramPrimitive *smp = new StringModelHistogramPrimitive ( theString );
  Material *histoMat = new Material;
  histoMat->ambient = vec4 ( 1,1,1,1 );
  histoMat->diffuse = vec4 ( 0,0,0,1 );
  histoMat->specular = vec4 ( 0,0,0,1 );
  histoMat->shininess = 0.0;
  histoMat->program = histoMat->loadShaders ( "MotionBlurredString" );
  smp->setMaterial(histoMat);
  // load the uniform array with the opacity function
  // XXX OMG this should be in a texture.
  float opacityValues[256];
  initOpacityValues ( opacityValues, 256 );
  glUseProgram ( histoMat->program );
  glUniform1fv(glGetUniformLocation(histoMat->program,"bins"), 256, opacityValues);
  glUseProgram(0);


  instance->addChild ( smp );


  FFTPrimitive *fft = new FFTPrimitive ( theString );
  Instance *fftTransform = new Instance;
  fftTransform->addChild ( fft );
  // on the guitar
  //  fftTransform->setMatrix ( glm::scale ( glm::translate ( glm::mat4(), 
  //							  glm::vec3 (-1.25, 0.1, 0.0) ),
  //					 glm::vec3(0.75, 1.0/128.0, 1.0) ) );

  // on the desk
  fftTransform->setMatrix ( glm::scale ( glm::translate ( glm::rotate ( glm::mat4(), -90.0f, glm::vec3(1.0,0.0,0.0) ),
							  glm::vec3 (-1.2, -0.3, -0.85) ),
					 glm::vec3(2.0, 1.0/128.0, 1.0) ) );
  instance->addChild ( fftTransform );

  Material *fftmat = new Material;
  //  fftmat->ambient = vec4 ( 1,1,1, 1.0 );
  fftmat->ambient = vec4 ( 0,0,0, 1.0 );
  fftmat->diffuse = vec4 ( 0,0,0, 1.0 );
  fftmat->specular = vec4 ( 0,0,0, 1.0 );
  fftmat->shininess = 0.0;
  fftmat->program = fftmat->loadShaders ( "PhongShading" );
  fftTransform->setMaterial(fftmat);


#ifdef USE_OSC
  // start the OSC thread
  if ( pthread_create ( &oscThreadId, 
			NULL, 
			oscThreadFunction,  
			static_cast<void *> (&oscParams) ) ) {
    cout << "error creating OSC thread!" << endl;
    exit(1);
  }
#endif
  
  //  std::cout << "MyWindow::init done" << std::endl;

  

}

void
MyWindow::oscCollectInput()
{
  map<std::string, float>::iterator i;
  for ( i = oscParams.value.begin(); i != oscParams.value.end(); i++ ) {
    guiSliders[i->first]->value ( i->second );
    guiSliders[i->first]->do_callback();
    oscParams.changed[i->first] = false;
  }
}

void 
MyWindow::draw() {
  if (!valid()) { 
    std::cout << "not valid!!!!" << std::endl;
    init(); 
    camera.setupPerspective( w(), h() );
  }
  root->update ( 0.033 );
  oscCollectInput();
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  camera.setDistance ( 1.0 );
  camera.draw(root);
}


int 
MyWindow::handle ( int event ) {
  //  std::cout << "handle " << event << ',' << Fl::event_x() << ',' << Fl::event_y() << std::endl;
  switch ( event ) {
  case FL_PUSH: // Fl::event_x() and Fl::event_y() 
    camera.startMouse ( Fl::event_x(), h()-Fl::event_y() );
    return 1;
  case FL_DRAG:
    camera.dragMouse ( Fl::event_x(), h()-Fl::event_y() );
    return 1;
  case FL_RELEASE:
    return 1;
  case FL_FOCUS:
  case FL_UNFOCUS:
    return 1; // to get kbd events
  case FL_KEYBOARD:
    // key in Fl::event(), ascii in Fl::event_text()
    // return 1 if understand/use the event
    //    std::cout << "FL_KEYBOARD " << Fl::event_key() << std::endl;
    keyboard ( Fl::event_key(), Fl::event_x(), Fl::event_y() );
    return 1;
  case FL_SHORTCUT:
    // key in Fl::event(), ascii in Fl::event_text()
    // return 1 if understand/use the event;
    return 1;
  default:
    return Fl_Gl_Window::handle(event);
  }
}

void 
idle (void *data) {
  Fl_Gl_Window *w = static_cast<Fl_Gl_Window*>(data);
  if (w)
    w->redraw();
}
///////////////////////////////////////////////////

int 
main(int argc, char** argv)
{
  // the enclosing FLTK window
  Fl_Window *window = new Fl_Window(winWidth, winHeight, "CSUEB Virtual String Apparatus");

  // the OpenGL subwindow
  MyWindow *mywindow = new MyWindow(0,0,winWidth,offsetWidgets,NULL);

  // make the gl window a sub-window
  window->add(mywindow);

  // gets some animation going
  Fl::add_idle(idle,mywindow);

  // add some controls
  Fl_Pack *bottomWindow = makeApparatusControls ( winWidth, winHeight, 
						  offsetWidgets + 20, 
						  coarsefineHeight );
  window->add ( bottomWindow );
  window->end();
  window->show(argc, argv);
  return Fl::run();
}


