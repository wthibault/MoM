#include "gui.h"
#include <algorithm>
#include "StringModel.h"
using namespace std;

extern StringModel *theString;

float 
getFloatValue ( void *input )
{
  return atof ( ((Fl_Float_Input*)input)->value () );
}

void 
setFloatValue ( void* input, float sliderValue ) 
{
  char str[100];
  sprintf(str, "%8.5f", sliderValue);
  ((Fl_Float_Input*)input)->value ( str );
}

void 
handleCoarse (Fl_Widget* o, void* input, float maxValue, float& lastFine) 
{
  Fl_Valuator *slider = ((Fl_Valuator*)o);
  float sliderValue = slider->value();
  float inputValue = getFloatValue ( input );
  float value;
  value = sliderValue * maxValue;
  value = max ( 0.0f, min ( maxValue, value ) );
  setFloatValue ( input, value );
}

void 
handleFine (Fl_Widget* o, void* input, float maxValue, float& lastFine) 
{
  Fl_Valuator *slider = ((Fl_Valuator*)o);
  float sliderValue = slider->value();
  float inputValue = getFloatValue ( input );
  float value;
  float incr = (sliderValue - lastFine) * 0.001;
  value = inputValue + incr;
  lastFine = sliderValue;
  value = max ( 0.0f, min ( maxValue, value ) );
  setFloatValue ( input, value );
}

//////////////////////
// callbacks
/////////////////////

static void sliderTensionCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 1.0, lastFine );
  float f = getFloatValue ( input );
  theString->Ktension = double(f);
}

static void sliderTensionFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 1.0, lastFine );
  float f = getFloatValue ( input );
  theString->Ktension = double(f);
}



static void sliderDampingCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 1.0, lastFine );
  float f = getFloatValue ( input );
  theString->Kdamping = double(f);

}

static void sliderDampingFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 1.0, lastFine );
  float f = getFloatValue ( input );
  theString->Kdamping = double(f);

}



static void sliderVibFreqCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 22100.0, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorFreq = double(f);

}

static void sliderVibFreqFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 22100.0, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorFreq = double(f);

}



static void sliderVibAmpCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 0.005, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorAmplitude = double(f);
}

static void sliderVibAmpFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 0.005, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorAmplitude = double(f);
}

static void inputTensionCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->Ktension = double(f);
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

static void vibOnCallback ( Fl_Widget* o )
{
  theString->vibratorPower ( ((Fl_Light_Button*) o)->value() );
}

static void sineButtonCallback ( Fl_Widget *o )
{
}

static void sawtoothButtonCallback ( Fl_Widget *o )
{
}


/////////////////////////////////////////////////////////////

Fl_Group *
makeCoarseFineControl ( int w, int h, const char *label, 
			Fl_Callback *fCoarseSlider, 
			Fl_Callback *fFineSlider,
			Fl_Callback *fInput )
{
  int inputHeight = 20;
  int inputWidth = 50;
  int horizMargin = 30;
  //  int horizOffset = w/5;
  int horizRemain = w - inputWidth - 2*horizMargin;
  int horizRemainEnd = horizRemain - horizMargin;

  Fl_Group *group = new Fl_Group ( 0,0, w, h, NULL );


  Fl_Float_Input *input = new Fl_Float_Input ( horizMargin,0, inputWidth, inputHeight, label );
  //  Fl_Float_Input *input = new Fl_Float_Input ( horizMargin,0, inputWidth, inputHeight, label );
  input->labelsize(10);
  input->callback ( fInput );
  group->add(input);

  Fl_Pack *sliderPack = new Fl_Pack ( inputWidth+horizMargin, 0, horizRemain,h );
  //  Fl_Hor_Nice_Slider *coarseSlider = new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, 2*h/3, "coarse" );
  Fl_Hor_Nice_Slider *coarseSlider = new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, 2*h/3 );
  coarseSlider->color ( FL_GRAY, FL_DARK_YELLOW );
  coarseSlider->callback ( fCoarseSlider, input );

  //  Fl_Hor_Nice_Slider *fineSlider =   new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, h/3, "fine" );
  Fl_Hor_Nice_Slider *fineSlider =   new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, h/3 );
  fineSlider->color ( FL_GRAY, FL_DARK_YELLOW );
  fineSlider->callback ( fFineSlider, input );

  sliderPack->add ( coarseSlider );
  sliderPack->add ( fineSlider );

  group->add(sliderPack);

  group->end();
  group->labelsize(10);
  return group;
}


Fl_Group* 
makeVibControls(int x, int y, int width, int height, int coarsefineHeight)
{
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  //  Fl_Pack *h = new Fl_Pack ( 0,0,winWidth, coarsefineHeight );
  //  h->type( Fl_Pack::HORIZONTAL );
  Fl_Group *h = new Fl_Group ( 0,0,width, coarsefineHeight );
  
  Fl_Light_Button *onbut = new Fl_Light_Button ( 0,0, 140, 20, "Vibrator on");
  onbut->callback ( vibOnCallback );
  h->add ( onbut );

  Fl_Button *o = new Fl_Light_Button (0,20, 140,20, "constant energy");
  o->tooltip("Turning on constant energy will cause the vibrator to have less amplitude at higher frequencies.  Turning it off makes it easier to find antinodes at higher frequencies.");
  h->add(o);

  Fl_Pack *buts = new Fl_Pack ( 150,0, 20, coarsefineHeight, "" );
  buts->type(Fl_Pack::VERTICAL);

  o = new Fl_Button(0,0, 20, 20, "sine");
  o->tooltip("Set vibrator waveform to sine");
  o->type(102);
  o->selection_color((Fl_Color)1);
  o->align(Fl_Align(FL_ALIGN_RIGHT));
  o->value(1);
  o->callback(sineButtonCallback);
  buts->add(o);

  o = new Fl_Button(0, 0, 20, 20, "sawtooth");
  o->tooltip("Set vibrator waveform to sawtooth");
  o->type(102);
  o->selection_color((Fl_Color)1);
  o->align(Fl_Align(FL_ALIGN_RIGHT));
  o->callback(sawtoothButtonCallback);
  buts->add(o);
  buts->end();

  h->add(buts);

  h->end();

  widgetPacker->add(h);

  Fl_Group *pack3 = makeCoarseFineControl(width,coarsefineHeight,"Vib. Freq.", 
					  sliderVibFreqCoarseCallback, 
					  sliderVibFreqFineCallback, 
					  inputVibFreqCallback);
  widgetPacker->add(pack3);
  Fl_Group *pack4 = makeCoarseFineControl(width,coarsefineHeight,"Vib. Amp.", 
					  sliderVibAmpCoarseCallback, 
					  sliderVibAmpFineCallback, 
					  inputVibAmpCallback );
  widgetPacker->add(pack4);
  widgetPacker->end();

  return widgetPacker;
}

Fl_Group* 
makeStringControls(int x, int y, int width, int height, int coarsefineHeight )
{
  //  Fl_Pack *widgetPacker = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  Fl_Group *pack1 = makeCoarseFineControl(width,coarsefineHeight,"Tension", 
					  sliderTensionCoarseCallback, 
					  sliderTensionFineCallback, 
					  inputTensionCallback );
  widgetPacker->add(pack1);
  Fl_Group *pack2 = makeCoarseFineControl(width,coarsefineHeight,"Damping", 
					  sliderDampingCoarseCallback, 
					  sliderDampingFineCallback, 
					  inputDampingCallback );
  widgetPacker->add(pack2);
  widgetPacker->end();

  return widgetPacker;
}

Fl_Pack* 
makeApparatusControls ( int winWidth, int winHeight, int offsetWidgets, int coarsefineHeight )
{
  Fl_Pack *packer = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  packer->add(makeStringControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets, coarsefineHeight));
  packer->add(makeVibControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets, coarsefineHeight));
  packer->add ( new Fl_Window ( winWidth/2, winHeight/2 ) ); // ?????
  packer->end();
  return packer;
}
