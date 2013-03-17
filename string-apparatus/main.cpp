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


using namespace glm;

ModelNode *root;          
Primitive *prim;          
//mat4 projectionMatrix;    
//mat4 modelviewMatrix;

int width, height; // window size (glut)

const int winWidth = 800;
const int winHeight = 600;
const int offsetWidgets = winHeight / 2;
const int coarsefineHeight = (winHeight - offsetWidgets) / 7;


Trackball trackball(320,240,240);


//const char *shader = "DepthMap";
//const char *shader = "ConstantShading";
//const char *shader = "BumpMappedTexturedPhongShading";
const char *shader = "PhongShading";

StringModel *theString;
RtAudio dac;
//UI ui("VeraSe.ttf");

///////////////////////////////////////////////////////////////////////


void 
reshape (int w, int h)
{
  width = w;
  height = h;
  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  //  projectionMatrix = perspective ( 65.0f, (GLfloat) w / (GLfloat) h, 1.0f, 100.0f );
  //  modelviewMatrix = translate ( mat4(), vec3(0,0,-2) );
  float halfw = float(w)/2.0f;
  float halfh = float(h)/2.0f;
  trackball = Trackball ( halfw, halfh, fmin ( halfw, halfw ) );
}

Instance *
createRandomInstanceNewSystem ()
{
  // create a new instance to refer to the same primitive, transformed
  ParticleSystem *prim = new ParticleSystem();
  
  Instance* anotherInstance = new Instance();
  anotherInstance->addChild ( prim );
  vec2 pos = vec2 ( urand() - 0.5f, urand() - 0.5f );
  float len = 2*glm::length(pos);
  anotherInstance->setMatrix ( rotate ( translate( mat4(), 2.0f * vec3 ( 0.0, pos.x, pos.y ) ),
					90.0f, vec3(0,0,1) ) );

  // create a material to use
  Material *mat = new Material;
  mat->ambient = vec4 ( 0.1,0.1,0.1, 1.0 );
  vec3 rimcolor = vec3 ( 1.0, 0.0, 0.0 );
  vec3 centercolor = vec3 ( 1.0, 1.0, 0.0 );
  vec3 color = (1-len)*centercolor + len*rimcolor;
  mat->diffuse = vec4 ( color, 1.0 );
  mat->specular = vec4 ( 1.0, 1.0, 1.0, 1.0 );
  mat->shininess = 300.0;
  mat->program = mat->loadShaders ( shader );
  anotherInstance->setMaterial ( mat );
  
  return anotherInstance;
}

void 
keyboard (unsigned char key, int x, int y)
{
  float upCoarse = 1.06;
  float upFine = 1.001;
  float downCoarse = 1 / upCoarse;
  float downFine = 1 / upFine;
  float dampAdjustment = 1e-5;

  switch (key) {
  case 'a':
    dynamic_cast<Instance*>(root)->addChild ( createRandomInstanceNewSystem() );
    break;

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
    theString->Kdamping += dampAdjustment;
    theString->Kdamping = fmin ( 1.0, theString->Kdamping );
    break;
  case 'E' : 
    // fine damping control - up
    //    theString->Kdamping *= upFine / dampAdjustment;
    theString->Kdamping += dampAdjustment/10;
    theString->Kdamping = fmin ( 1.0, theString->Kdamping );
    break;
  case 'd' : 
    // coarse damping control - down
    //    theString->Kdamping *= downCoarse / dampAdjustment;
    theString->Kdamping -= dampAdjustment;
    break;
  case 'D' : 
    // fine damping control - down
    //    theString->Kdamping *= downFine / dampAdjustment;
    theString->Kdamping -= dampAdjustment/10;
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
    //    values[i] = min ( 1.0, 1.0 / fmax ( EPSILON, M_PI * sqrt ( 1 - (2*x-1) * (2*x-1) ) ) );
    values[i] = min ( 1.0, 1.0 / (M_PI * sqrt ( 1 - (2*x-1) * (2*x-1) )) );
    //    std::cout << values[i] << std::endl;
  }
}

void 
init (int argc, char **argv)
{
  
  // XXX this should pull it apart so we can move things
  if ( argc == 2 ) {
    int choice = atoi ( argv[1] );
    switch (choice) {
    case 1:  
      prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );
      break;
    case 2:
      prim = new ObjFilePrimitive ( "objfiles/string-scene.obj" );
      break;
    default:
      prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );
      break;
    }
  } else {
      prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );
  }

  // create a root Instance to contain this primitive
  Instance *instance = new Instance();
  instance->setMatrix ( mat4() );
  instance->addChild ( prim );

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

  if (argc >= 2 )
    mat->program = mat->loadShaders ( "PhongShading" );
  else
    mat->program = mat->loadShaders ( shader );

  // attach the material to the background object
  //  instance->setMaterial ( mat );

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

  // audio params
  int sampleRate = 44100;
  unsigned int bufferFrames = 256; // 256 sample frames ~ 5ms 
  //unsigned int bufferFrames = 1024; // 256 sample frames ~ 5ms 

  // the simulation 
  theString = new StringModel ( NUM_MASSES, 0.01, 0.99999, 2, sampleRate, bufferFrames );

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
#ifdef USE_HISTOGRAMS
  StringModelHistogramPrimitive *smp = new StringModelHistogramPrimitive ( theString );
  Material *histoMat = new Material;
  histoMat->ambient = vec4 ( 1,1,1,1 );
  histoMat->diffuse = vec4 ( 0,0,0,1 );
  histoMat->specular = vec4 ( 0,0,0,1 );
  histoMat->shininess = 0.0;
  histoMat->program = histoMat->loadShaders ( "MotionBlurredString" );
  smp->setMaterial(histoMat);
  // load the uniform array with the opacity function
  // OMG this should be in a texture.
  float opacityValues[256];
  initOpacityValues ( opacityValues, 256 );
  glUseProgram ( histoMat->program );
  glUniform1fv(glGetUniformLocation(histoMat->program,"bins"), 256, opacityValues);
  glUseProgram(0);

#else
  StringModelPrimitive *smp = new StringModelPrimitive ( theString );
#endif

  instance->addChild ( smp );

  FFTPrimitive *fft = new FFTPrimitive ( theString );
  Instance *fftTransform = new Instance;
  fftTransform->addChild ( fft );
  fftTransform->setMatrix ( glm::scale ( glm::translate ( glm::mat4(), 
							  glm::vec3 (-1.25, 0.1, 0.0) ),
					 glm::vec3(0.75, 1.0/128.0, 1.0) ) );
  instance->addChild ( fftTransform );

  Material *fftmat = new Material;
  //  fftmat->ambient = vec4 ( 1,1,1, 1.0 );
  fftmat->ambient = vec4 ( 0,0,0, 1.0 );
  fftmat->diffuse = vec4 ( 0,0,0, 1.0 );
  fftmat->specular = vec4 ( 0,0,0, 1.0 );
  fftmat->shininess = 0.0;
  fftmat->program = fftmat->loadShaders ( "PhongShading" );
  fftTransform->setMaterial(fftmat);

}

void 
mouse ( int button, int state, int x, int y )
{
  if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {
    trackball.startMouse ( x, height-y );
  }
}

void
motion ( int x, int y ) 
{
  trackball.dragMouse(x,height-y);
}
/////////////////////////////////////////////////////
class MyWindow : public Fl_Gl_Window {
  void draw();
  int handle(int);
  ModelNode *root;
  Camera     camera;
public:
  MyWindow ( int x, int y , int w, int h, const char *L );
  void init();
};

MyWindow::MyWindow ( int x, int y , int w, int h, const char *L )
    : Fl_Gl_Window ( x, y, w, h, L ) 
{
  mode ( FL_RGB8 | FL_DOUBLE | FL_ALPHA | FL_DEPTH | FL_MULTISAMPLE );
  end();
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


  Primitive *prim;
  prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );

  // create a root Instance to contain this primitive
  Instance *instance = new Instance();
  instance->setMatrix ( mat4() );
  instance->addChild ( prim );

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

  //  if (argc >= 2 )
  //    mat->program = mat->loadShaders ( "PhongShading" );
  //  else
  mat->program = mat->loadShaders ( shader );

  // attach the material to the background object
  //  instance->setMaterial ( mat );

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

  // audio params
  int sampleRate = 44100;
  unsigned int bufferFrames = 256; // 256 sample frames ~ 5ms 
  //unsigned int bufferFrames = 1024; // 256 sample frames ~ 5ms 

  // the simulation 
  theString = new StringModel ( NUM_MASSES, 0.01, 0.99999, 2, sampleRate, bufferFrames );

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
#ifdef USE_HISTOGRAMS
  StringModelHistogramPrimitive *smp = new StringModelHistogramPrimitive ( theString );
  Material *histoMat = new Material;
  histoMat->ambient = vec4 ( 1,1,1,1 );
  histoMat->diffuse = vec4 ( 0,0,0,1 );
  histoMat->specular = vec4 ( 0,0,0,1 );
  histoMat->shininess = 0.0;
  histoMat->program = histoMat->loadShaders ( "MotionBlurredString" );
  smp->setMaterial(histoMat);
  // load the uniform array with the opacity function
  // OMG this should be in a texture.
  float opacityValues[256];
  initOpacityValues ( opacityValues, 256 );
  glUseProgram ( histoMat->program );
  glUniform1fv(glGetUniformLocation(histoMat->program,"bins"), 256, opacityValues);
  glUseProgram(0);

#else
  StringModelPrimitive *smp = new StringModelPrimitive ( theString );
#endif

  instance->addChild ( smp );

  FFTPrimitive *fft = new FFTPrimitive ( theString );
  Instance *fftTransform = new Instance;
  fftTransform->addChild ( fft );
  fftTransform->setMatrix ( glm::scale ( glm::translate ( glm::mat4(), 
							  glm::vec3 (-1.25, 0.1, 0.0) ),
					 glm::vec3(0.75, 1.0/128.0, 1.0) ) );
  instance->addChild ( fftTransform );

  Material *fftmat = new Material;
  //  fftmat->ambient = vec4 ( 1,1,1, 1.0 );
  fftmat->ambient = vec4 ( 0,0,0, 1.0 );
  fftmat->diffuse = vec4 ( 0,0,0, 1.0 );
  fftmat->specular = vec4 ( 0,0,0, 1.0 );
  fftmat->shininess = 0.0;
  fftmat->program = fftmat->loadShaders ( "PhongShading" );
  fftTransform->setMaterial(fftmat);

}

void MyWindow::draw() {
  if (!valid()) { 
    std::cout << "not valid!!!!" << std::endl;
    init(); // XXX not on reshape, just once please. in constructor?
    camera.setupPerspective( w(), h() );
  }

  // draw

  camera.draw(root);
  draw_children();
}

int MyWindow::handle ( int event ) {
  //  std::cout << "handle " << event << std::endl;
  switch ( event ) {
  case FL_PUSH: // Fl::event_x() and Fl::event_y() 
    return 1;
  case FL_DRAG:
    return 1;
  case FL_RELEASE:
    return 1;
  case FL_FOCUS:
  case FL_UNFOCUS:
    return 1; // to get kbd events
  case FL_KEYBOARD:
    // key in Fl::event(), ascii in Fl::event_text()
    // return 1 if understand/use the event;
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
  Fl_Widget *w = static_cast<Fl_Widget*>(data);
  if (w)
    w->redraw();
}
///////////////////////////////////////////////////

int 
main(int argc, char** argv)
{

  // the enclosing FLTK window
  Fl_Window *window = new Fl_Window(winWidth, winHeight);

  // the OpenGL subwindow
  MyWindow *mywindow = new MyWindow(0,0,winWidth,offsetWidgets,NULL);

  // make the gl window a sub-window
  window->add(mywindow);

  // gets some animation going
  Fl::add_idle(idle,mywindow);

  // add some controls
  Fl_Pack *bottomWindow = makeApparatusControls ( winWidth, winHeight, offsetWidgets, coarsefineHeight );
  window->add ( bottomWindow );

  window->end();

  window->show(argc, argv);

  return Fl::run();

  

  return 0;
}


