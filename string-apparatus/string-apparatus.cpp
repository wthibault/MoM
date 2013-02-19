//
// example3.cpp
// particle system, trackball
//

#include "ssg.h"
#include <iomanip>
#include "ParticleSystem.h"
#include "Trackball.h"
#include "StringModel.h"
#include <fftw3.h>
#include "ui.h"

using namespace glm;

ModelNode *root;          
Primitive *prim;          
mat4 projectionMatrix;    
mat4 modelviewMatrix;

int width, height; // window size

Trackball trackball(320,240,240);


//const char *shader = "DepthMap";
//const char *shader = "ConstantShading";
const char *shader = "BumpMappedTexturedPhongShading";

StringModel *theString;
RtAudio dac;
UI ui("VeraSe.ttf");

//////
/////////////////////////////////////////////////////////////////////////
// StringModelPrimitive - an SSG rendering Primitive subclassed from ParticleSystem
//

class StringModelPrimitive : public ParticleSystem
{
public:
  StringModelPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel ), renderScale_ ( 10.0 )
  {
    StringModelPrimitive::init();
  }
  ~StringModelPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      points_.push_back ( glm::vec3(0,0,0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( i );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 0) );
    }
    drawingPrimitive_ = GL_POINTS;
    Primitive::init();
  }

  void update(float dt) 
  {
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // copy the positions of masses to the vertex array
    glBindVertexArray ( vao_ );

    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    float deltaX = 1.0 / theString_->numMasses;
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      float x,y,z;
      x = 2*(i * deltaX)-1;
      y = renderScale_ * theString_->yold[i];
      z = 0;
      points_.push_back( glm::vec3 ( x,y,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      indices_.push_back ( i );
    }
    
    long int sizeofPoints = sizeof(glm::vec3)*points_.size();
    int sizeofNormals = sizeof(glm::vec3)*normals_.size();
    int sizeofTexCoords = sizeof(glm::vec2)*texCoords_.size();
    drawingPrimitive_ = GL_POINTS;

    glBindBuffer ( GL_ARRAY_BUFFER, arrayBuffer_ );
    glBufferSubData( GL_ARRAY_BUFFER, 0, sizeofPoints, &points_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints, sizeofNormals, &normals_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints + sizeofNormals, sizeofTexCoords, &texCoords_[0] );
    glBindBuffer ( GL_ARRAY_BUFFER, 0 );
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer_);
    int sizeofIndices = indices_.size()*sizeof(unsigned int);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeofIndices, &indices_[0]);
    glBindBuffer ( GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    
    Primitive::draw ( mv, proj, mat );

    //    printParams();
  }

  void printParam ( const std::string name, float value ) {
    std::cout << name << std::setw( 6 ) << value << std::endl;
  }

  void printParams() {
    // temp until HUD
    std::cout << std::endl << std::endl << std::endl;
    printParam ("Ktension    ", theString_->Ktension );
    printParam ("Kdamping    ", theString_->Kdamping );
    printParam ("vib freq.   ", theString_->vibratorFreq );
    printParam ("vib amp.    ", theString_->vibratorAmplitude );
    
    std::cout << std::endl;
    std::cout << "t/g tension (coarse) " << std::endl;
    std::cout << "T/G tension (fine) " << std::endl;
    std::cout << "e/d damping (coarse) " << std::endl;
    std::cout << "E/D damping (fine) " << std::endl;
    std::cout << "u/j vib freq (coarse) " << std::endl;
    std::cout << "U/J vib freq (fine) " << std::endl;
    std::cout << "i/k vib amp (coarse) " << std::endl;
    std::cout << "I/K vib amp (fine) " << std::endl;
    std::cout << "p/P pluck" << std::endl;
    std::cout << "v toggle vibrator " << std::endl;
    std::cout << "r reset" << std::endl;
    std::cout << "ESC exit" << std::endl;
  }

  StringModel *theString_;
  float        renderScale_;
};


/////////////////////////////////////////////////////////////////

class StringModelHistogramPrimitive : public ParticleSystem
{
public:
  StringModelHistogramPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel ), renderScale_ ( 10.0 )
  {
    StringModelHistogramPrimitive::init();
  }
  ~StringModelHistogramPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      points_.push_back ( glm::vec3(0,0,0) );
      points_.push_back ( glm::vec3(0,0,0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( i );
      indices_.push_back ( i );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 0) );
      texCoords_.push_back ( glm::vec2(i/theString_->numMasses, 1) );
    }
    drawingPrimitive_ = GL_LINES;
    Primitive::init();
  }

  void update(float dt) 
  {
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // draw a line based on histogram of positions
    glBindVertexArray ( vao_ );

    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    float deltaX = 1.0 / theString_->numMasses;
    for ( int i = 0; i < theString_->numMasses; i++ ) {
      float x,y0,y1,z;
      x = 2*(i * deltaX)-1;
      y0 = renderScale_ * theString_->histograms[i].minVal;
      y1 = renderScale_ * theString_->histograms[i].maxVal;
      z = 0;
      points_.push_back( glm::vec3 ( x,y0,z ) );
      points_.push_back( glm::vec3 ( x,y1,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
      theString_->histograms[i].clear();
      theString_->histograms[i].minVal = -1e-6;
      theString_->histograms[i].maxVal = +1e-6;
    }

    
    long int sizeofPoints = sizeof(glm::vec3)*points_.size();
    int sizeofNormals = sizeof(glm::vec3)*normals_.size();
    int sizeofTexCoords = sizeof(glm::vec2)*texCoords_.size();
    drawingPrimitive_ = GL_LINES;

    glBindBuffer ( GL_ARRAY_BUFFER, arrayBuffer_ );
    glBufferSubData( GL_ARRAY_BUFFER, 0, sizeofPoints, &points_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints, sizeofNormals, &normals_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints + sizeofNormals, sizeofTexCoords, &texCoords_[0] );
    glBindBuffer ( GL_ARRAY_BUFFER, 0 );
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer_);
    int sizeofIndices = indices_.size()*sizeof(unsigned int);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeofIndices, &indices_[0]);
    glBindBuffer ( GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    
    Primitive::draw ( mv, proj, mat );

    //    printParams();
  }

  void printParam ( const std::string name, float value ) {
    std::cout << name << std::setw( 6 ) << value << std::endl;
  }

  void printParams() {
    // temp until HUD
    std::cout << std::endl << std::endl << std::endl;
    printParam ("Ktension    ", theString_->Ktension );
    printParam ("Kdamping    ", theString_->Kdamping );
    printParam ("vib freq.   ", theString_->vibratorFreq );
    printParam ("vib amp.    ", theString_->vibratorAmplitude );
    
    std::cout << std::endl;
    std::cout << "t/g tension (coarse) " << std::endl;
    std::cout << "T/G tension (fine) " << std::endl;
    std::cout << "e/d damping (coarse) " << std::endl;
    std::cout << "E/D damping (fine) " << std::endl;
    std::cout << "u/j vib freq (coarse) " << std::endl;
    std::cout << "U/J vib freq (fine) " << std::endl;
    std::cout << "i/k vib amp (coarse) " << std::endl;
    std::cout << "I/K vib amp (fine) " << std::endl;
    std::cout << "p/P pluck" << std::endl;
    std::cout << "v toggle vibrator " << std::endl;
    std::cout << "r reset" << std::endl;
    std::cout << "ESC exit" << std::endl;
  }

  StringModel *theString_;
  float        renderScale_;
};





/////////////////////////////////////////////////////////////////////////


class FFTPrimitive : public ParticleSystem
{
public:
  FFTPrimitive( StringModel *stringModel ) 
    : theString_ ( stringModel )
  {
    FFTPrimitive::init();
  }
  ~FFTPrimitive () {}

  void init() 
  {
    points_.clear();
    normals_.clear();
    indices_.clear();
    texCoords_.clear();
    int nFreqBins = theString_->numFramesToAnalyze / 2;
    float step = 1.0 / nFreqBins;
    for ( int i = 0; i < nFreqBins; i++ ) {
      points_.push_back ( glm::vec3(i*step, 0.0, 0) );
      points_.push_back ( glm::vec3(i*step, 1.0, 0) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
      texCoords_.push_back ( glm::vec2(i*step, 0) );
      texCoords_.push_back ( glm::vec2(i*step, 1) );
    }
    drawingPrimitive_ = GL_LINES;
    Primitive::init();
  }

  void update(float dt) 
  {
    // setup the fft
    // XXX move out of string into here XXX
    // copy them backwards, same diff to the fft.
    theString_->ringBuffer.copyReversed ( theString_->fftwIn, theString_->numFramesToAnalyze, 2 );

    // run the fft on them
    fftw_execute ( theString_->fftwPlan );

    // load the geometry into the vao
    points_.clear();
    indices_.clear();
    normals_.clear();
    texCoords_.clear();

    int nFreqBins = theString_->numFramesToAnalyze / 2;
    float deltaX = 1.0 / nFreqBins;
    for ( int i = 0; i < nFreqBins; i++ ) {
      float x,y,z;
      x = i * deltaX;
      y = 20*log10(fabs(theString_->fftwOut[i][0]));
      z = 0;
      points_.push_back( glm::vec3 ( x,-100,z ) );
      points_.push_back( glm::vec3 ( x,y,z ) );
      normals_.push_back ( glm::vec3(0,0,1) );
      normals_.push_back ( glm::vec3(0,0,1) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 0) );
      texCoords_.push_back ( glm::vec2(i * deltaX, 1) );
      indices_.push_back ( 2*i );
      indices_.push_back ( 2*i+1 );
    }
  }

  void draw ( glm::mat4 mv, glm::mat4 proj, Material *mat ) 
  {
    // copy the positions of masses to the vertex array
    glBindVertexArray ( vao_ );

    long int sizeofPoints = sizeof(glm::vec3)*points_.size();
    int sizeofNormals = sizeof(glm::vec3)*normals_.size();
    int sizeofTexCoords = sizeof(glm::vec2)*texCoords_.size();
    drawingPrimitive_ = GL_LINES;

    glBindBuffer ( GL_ARRAY_BUFFER, arrayBuffer_ );
    glBufferSubData( GL_ARRAY_BUFFER, 0, sizeofPoints, &points_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints, sizeofNormals, &normals_[0] );
    glBufferSubData( GL_ARRAY_BUFFER, sizeofPoints + sizeofNormals, sizeofTexCoords, &texCoords_[0] );
    glBindBuffer ( GL_ARRAY_BUFFER, 0 );
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementBuffer_);
    int sizeofIndices = indices_.size()*sizeof(unsigned int);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeofIndices, &indices_[0]);
    glBindBuffer ( GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glBindVertexArray(0);
    
    Primitive::draw ( mv, proj, mat );

    //    printParams();
  }

  StringModel *theString_;
  
};
///////////////////////////////////////////////////////////////////////

void 
display ()
{
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  root->update(0.033); // assumes called at 30Hz

  // mouse control 
  mat4 tball = trackball.getMat4();
  mat4 mv = modelviewMatrix * tball;

  root->draw(mv, projectionMatrix );

  ui.draw( width, height, theString );

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
init (int argc, char **argv)
{
  
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
  options.flags |= RTAUDIO_HOG_DEVICE;

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
#else
  StringModelPrimitive *smp = new StringModelPrimitive ( theString );
#endif

  instance->addChild ( smp );

  FFTPrimitive *fft = new FFTPrimitive ( theString );
  Instance *fftTransform = new Instance;
  fftTransform->addChild ( fft );
  fftTransform->setMatrix ( glm::scale ( glm::translate ( glm::mat4(), glm::vec3 (-1.5,-1.5,0.5) ),
					 glm::vec3(2,1/128.0,1) ) );
  instance->addChild ( fftTransform );

  Material *fftmat = new Material;
  fftmat->ambient = vec4 ( 0.1, 0.1, 0.2, 1.0 );
  fftmat->diffuse = vec4 ( 0.3, 0.3, 1.0, 1.0 );
  fftmat->specular = vec4 ( 1.0, 1.0, 1.0, 1.0 );
  fftmat->shininess = 133.0;
  fft->setMaterial(fftmat);

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
  glutTimerFunc(REFRESH_PERIOD,timer,REFRESH_PERIOD); 
  glutMainLoop();
  return 0;
}


