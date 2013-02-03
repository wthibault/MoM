//
// example3.cpp
// particle system, trackball
//

#include "ssg.h"
#include "ParticleSystem.h"
#include "Trackball.h"
#include "StringModel.h"
using namespace glm;

ModelNode *root;          
Primitive *prim;          
mat4 projectionMatrix;    
mat4 modelviewMatrix;

int width, height; // window size

Trackball trackball(320,240,240);


//const char *shader = "DepthMap";
const char *shader = "ConstantShading";

StringModel *theString;
RtAudio dac;

void 
display ()
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  root->update(0.033); // assumes called at 30Hz

  // mouse control 
  mat4 tball = trackball.getMat4();
  mat4 mv = modelviewMatrix * tball;

  root->draw(mv, projectionMatrix );

  glutSwapBuffers();
}

void 
timer ( int delay )
{
  glutTimerFunc( delay, timer, delay );
  glutPostRedisplay();
}

void 
reshape (int w, int h)
{
  width = w;
  height = h;
  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  projectionMatrix = perspective ( 65.0f, (GLfloat) w / (GLfloat) h, 1.0f, 100.0f );
  modelviewMatrix = translate ( mat4(), vec3(0,0,-4) );
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
  switch (key) {
  case 'a':
    dynamic_cast<Instance*>(root)->addChild ( createRandomInstanceNewSystem() );
    break;

  case 't' : theString->Ktension *= 1.05946; 
    theString->Ktension = fmin ( 1.0, theString->Ktension );
    break;
  case 'l' : theString->Ktension *= 0.943876; break;
  case 'P' : theString->pluck(); break;
  case 'p' : theString->pluckvel(); break;
  case 'r' : theString->reset(); break;
  case 'd' : theString->print(); break;
  case 'v' : theString->toggleVibrator(); break;
  case 'f' : theString->vibratorFreq *= 0.943875; 
    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;
  case 'F' : theString->vibratorFreq /= 0.943875; 
    std::cout << "vib freq = " << theString->vibratorFreq << std::endl;
    break;

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
init (int argc, char **argv)
{
  
  // create a primitive.  if supplied on command line, read a .obj wavefront file

  // if ( argc >= 2 ) {
  //   prim = new ObjFilePrimitive ( argv[1] );
  // } else {
  //   //    prim = new Triangle;
  //   prim = new ParticleSystem();
  //   std::cout << "usage: " << argv[0] << " [objfile.obj]\n";
  // }
  // std::cout << "Hit 'a' to add a new particle system.\n";

  // XXX this should pull it apart so we can move things
  prim = new ObjFilePrimitive ( "objfiles/string-apparatus.obj" );
  // create a root Instance to contain this primitive
  Instance *instance = new Instance();
  instance->setMatrix ( mat4() );
  instance->addChild ( prim );

  // the lights are global for all objects in the scene
  RenderingEnvironment::getInstance().lightPosition = vec4 ( 0,0,10,1 );
  RenderingEnvironment::getInstance().lightColor = vec4 ( 1,1,1,1 );
  // gravity center is in world coords, used by all ParticleSystem instances
  ParticleSystem::gravityCenter = vec3 ( 0,-1.9,0 );

  // create a material to use
  Material *mat = new Material;
  mat->ambient = vec4 ( 0.1, 0.1, 0.2, 1.0 );
  mat->diffuse = vec4 ( 0.5, 0.5, 0.1, 1.0 );
  mat->specular = vec4 ( 1.0, 1.0, 1.0, 1.0 );
  mat->shininess = 133.0;

  if (argc >= 2 )
    mat->program = mat->loadShaders ( "PhongShading" );
  else
    mat->program = mat->loadShaders ( shader );

  // attach the material to the primitive
  instance->setMaterial ( mat );

  // set the instance as the scene root
  root = instance;

  // misc OpenGL state
  glClearColor (0.0, 0.0, 0.0, 1.0);
  glEnable(GL_DEPTH_TEST);
  glPointSize(1.0);
  glLineWidth(1.0);

  // the simulation 
  theString = new StringModel ( 1000, 0.5, 0.99999, 8 );

  // the audio
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
    
  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  int sampleRate = 44100;
  unsigned int bufferFrames = 256; // 256 sample frames ~ 5ms 

  try { 
    dac.openStream ( &parameters, 
		     NULL, 
		     RTAUDIO_FLOAT32,
		     sampleRate, 
		     &bufferFrames, 
		     StringModel::audioCallback,
		     (void *)theString );
    dac.startStream();
  }
  catch ( RtError& e ) {
    std::cout << "\nexception on dac:\n";
    e.printMessage();
    exit(0);
  }
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


int 
main(int argc, char** argv)
{
  glutInit(&argc, argv);
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize (300, 300); 
  glutInitWindowPosition (100, 100);
  glutCreateWindow (argv[0]);
  
#ifndef __APPLE__
  glewInit();
#endif
  init(argc,argv);
  
  glutMouseFunc ( mouse );
  glutMotionFunc ( motion );
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutTimerFunc(30,timer,30); 
  glutMainLoop();
  return 0;
}


