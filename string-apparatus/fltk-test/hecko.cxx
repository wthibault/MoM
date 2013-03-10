#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Hor_Nice_Slider.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Toggle_Button.H>
#include <iostream>
#ifndef __APPLE__
#include <GL/glew.h>
#endif
#include "ssg.h"
#include "Camera.h"
using namespace glm;

const int winWidth = 800;
const int winHeight = 600;
const int offsetWidgets = winHeight / 2;
const int coarsefineHeight = (winHeight - offsetWidgets) / 7;

class MyWindow : public Fl_Gl_Window {
  void draw();
  int handle(int);
  ModelNode *root;
  Camera     camera;
public:
  MyWindow ( int x, int y , int w, int h, const char *L );
  //    : Fl_Gl_Window ( x, y, w, h, L ) {
  void init();
};

MyWindow::MyWindow ( int x, int y , int w, int h, const char *L )
    : Fl_Gl_Window ( x, y, w, h, L ) 
{
  end();
}

void
MyWindow::init()
{

  CheckError();
#ifndef __APPLE__
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if ( err != GLEW_OK ) {
      std::cout << "glewInit() fail, bailing." << std::endl;
      exit(1);
    }
#endif

    //  create a primitive.  if supplied on command line, read a .obj wavefront file
  Primitive *prim = new Triangle;

  // create the graph
  // the second child of the root must be an instance: it will have its matrix changed
  Instance *scene = new Instance();
  scene->setMatrix ( mat4() );
  scene->addChild ( prim );
  Instance *mover = new Instance();
  mover->addChild ( prim );
  scene->addChild ( mover );

  // enable camera trackball
  camera.enableTrackball (true);

  // the lights are global for all objects in the scene
  RenderingEnvironment::getInstance().lightPosition = vec4 ( 4,10,5,1 );
  RenderingEnvironment::getInstance().lightColor = vec4 ( 1,1,1,1 );


  // create a material to use
  Material *mat = new Material;
  mat->ambient = vec4 ( 0.1, 0.1, 0.2, 1.0 );
  mat->diffuse = vec4 ( 0.5, 0.5, 0.1, 1.0 );
  mat->specular = vec4 ( 1.0, 1.0, 1.0, 1.0 );
  mat->shininess = 133.0;
  mat->program = mat->loadShaders ( "PhongShading" );

  // attach the material to the primitive
  scene->setMaterial ( mat );

  // set the instance as the scene root
  root = scene;

  // setup projection using w() and h()
  glClearColor (0.0, 0.0, 0.0, 1.0);
  glEnable(GL_DEPTH_TEST);
  CheckError();
}

void MyWindow::draw() {
  if (!valid()) { 
    init();
    camera.setupPerspective( w(), h() );
  }
  // draw
  //  std::cout << "draw" << std::endl;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  camera.draw(root);
  draw_children();
}

int MyWindow::handle ( int event ) {
  std::cout << "handle " << event << std::endl;
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
void idle (void *data) {
  //  std::cout << "idle\n";
  Fl_Widget *w = static_cast<Fl_Widget*>(data);
  if (w)
    w->redraw();
}

float 
getFloatValue ( void *input )
{
  return atof ( ((Fl_Float_Input*)input)->value () );
}

void 
setFloatValue ( void* input, float sliderValue ) 
{
  char str[100];
  sprintf(str, "%f", sliderValue);
  ((Fl_Float_Input*)input)->value ( str );
}

void 
handleCoarseFine (Fl_Widget* o, void* input, float maxValue, float& lastFine) 
{
  Fl_Valuator *slider = ((Fl_Valuator*)o);
  float sliderValue = slider->value();
  float inputValue = getFloatValue ( input );
  float value;
  if ( strcmp ( slider->label(), "coarse" ) == 0 ) {
    value = sliderValue * maxValue;
  } else {
    float incr = (sliderValue - lastFine) * 0.001;
    value = inputValue + incr;
    lastFine = sliderValue;
  }

  value = max ( 0.0f, min ( maxValue, value ) );
  setFloatValue ( input, value );
}

static void sliderTensionCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarseFine ( o, input, 1.0, lastFine );
}

static void sliderDampingCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarseFine ( o, input, 1.0, lastFine );

}

static void sliderVibFreqCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarseFine ( o, input, 22100.0, lastFine );
}

static void sliderVibAmpCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarseFine ( o, input, 5.0, lastFine );
}

static void inputTensionCallback ( Fl_Widget* o, void *theFloat )
{
}

static void inputDampingCallback ( Fl_Widget* o, void *theFloat )
{
}

static void inputVibFreqCallback ( Fl_Widget* o, void *theFloat )
{
}

static void inputVibAmpCallback ( Fl_Widget* o, void *theFloat )
{
}


Fl_Group *
makeCoarseFineControl ( int w, int h, const char *label, Fl_Callback *fSlider, Fl_Callback *fInput )
{
  int horizMargin = 50;
  int horizOffset = w/5;
  int horizRemain = w - horizOffset;
  int horizRemainEnd = w - horizMargin;
  int inputHeight = 20;

  Fl_Group *group = new Fl_Group ( 0,0, w, h, NULL );


  Fl_Float_Input *input = new Fl_Float_Input ( horizMargin,0, horizOffset-horizMargin, inputHeight, label );
  input->labelsize(10);
  input->callback ( fInput );
  group->add(input);

  Fl_Pack *sliderPack = new Fl_Pack ( horizOffset, 0, horizRemain,h );
  Fl_Hor_Nice_Slider *coarseSlider = new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, 2*h/3, "coarse" );
  coarseSlider->color ( FL_GRAY, FL_DARK_YELLOW );
  coarseSlider->callback ( fSlider, input );

  Fl_Hor_Nice_Slider *fineSlider =   new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, h/3, "fine" );
  fineSlider->color ( FL_GRAY, FL_DARK_YELLOW );
  fineSlider->callback ( fSlider, input );

  sliderPack->add ( coarseSlider );
  sliderPack->add ( fineSlider );

  group->add(sliderPack);

  group->end();
  group->labelsize(10);
  return group;
}


Fl_Group* makeVibControls(int x, int y, int width, int height)
{
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  Fl_Pack *h = new Fl_Pack ( 0,0,winWidth, coarsefineHeight );
  h->type( Fl_Pack::HORIZONTAL );
  
  Fl_Light_Button *onbut = new Fl_Light_Button ( 0,0, 100, 20, "Vibrator on" );
  h->add ( onbut );

  Fl_Pack *buts = new Fl_Pack ( 0,0, 20, 100, "waveform" );
  buts->type(Fl_Pack::VERTICAL);

  Fl_Button* o = new Fl_Button(0,0, 20, 20, "sine");
  o->tooltip("Set vibrator waveform to sine");
  o->type(102);
  o->selection_color((Fl_Color)1);
  o->align(Fl_Align(FL_ALIGN_RIGHT));
  buts->add(o);

  o = new Fl_Button(0,0, 20, 20, "sawtooth");
  o->tooltip("Set vibrator waveform to sawtooth");
  o->type(102);
  o->selection_color((Fl_Color)1);
  o->align(Fl_Align(FL_ALIGN_RIGHT));
  buts->add(o);
  buts->end();

  h->add(buts);

  o = new Fl_Toggle_Button (0,0, 100,20, "constant energy");
  o->tooltip("Turning on constant energy will cause the vibrator to have less amplitude at higher frequencies.  Turning it off makes it easier to find antinodes at higher frequencies.");
  h->add(o);

  widgetPacker->add(h);



  Fl_Group *pack3 = makeCoarseFineControl(winWidth,coarsefineHeight,"Vib. Freq.", 
					  sliderVibFreqCallback, inputVibFreqCallback);
  widgetPacker->add(pack3);
  Fl_Group *pack4 = makeCoarseFineControl(winWidth,coarsefineHeight,"Vib. Amp.", 
					  sliderVibAmpCallback, inputVibAmpCallback );
  widgetPacker->add(pack4);
  widgetPacker->end();

  return widgetPacker;
}

Fl_Group* makeStringControls(int x, int y, int width, int height)
{
  //  Fl_Pack *widgetPacker = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  Fl_Group *pack1 = makeCoarseFineControl(winWidth,coarsefineHeight,"Tension", 
					  sliderTensionCallback, inputTensionCallback );
  widgetPacker->add(pack1);
  Fl_Group *pack2 = makeCoarseFineControl(winWidth,coarsefineHeight,"Damping", 
					  sliderDampingCallback, inputDampingCallback );
  widgetPacker->add(pack2);
  widgetPacker->end();

  return widgetPacker;
}


int 
main(int argc, char **argv) {

  // the enclosing FLTK window
  Fl_Window *window = new Fl_Window(winWidth, winHeight);

  // the OpenGL subwindow
  MyWindow *mywindow = new MyWindow(0,0,winWidth,offsetWidgets,NULL);

  // make the gl window a sub-window
  window->add(mywindow);

  // gets some animation going
  Fl::add_idle(idle,mywindow);

  // pack some controls

  Fl_Pack *bottomWindow = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  bottomWindow->add(makeStringControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets));
  bottomWindow->add(makeVibControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets));
  bottomWindow->add ( new Fl_Window ( winWidth/2, winHeight/2 ) ); // ?????
  bottomWindow->end();

 
 window->add ( bottomWindow );
  window->end();
  window->show(argc, argv);
  return Fl::run();
}
