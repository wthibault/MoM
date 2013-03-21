#include "gui.h"
#include <algorithm>
#include "StringModel.h"
using namespace std;

extern StringModel *theString;
extern const double initVibFreq;
extern const double initVibAmp;

#ifdef NEW_STRING_MODEL
extern const double initHangerMass;
extern const double initDecayTime;
extern const double initMassDensity;
#else
extern const double initTension;
extern const double initDamping;
#endif

const float maxVibAmp = 0.005;
const float incrFine = 0.1;

const char *floatFormat = "%8.5f";

float 
getFloatValue ( void *input )
{
  return atof ( ((Fl_Float_Input*)input)->value () );
}

void 
setFloatValue ( void* input, float sliderValue ) 
{
  char str[100];
  sprintf(str, floatFormat, sliderValue);
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
  float incr = (sliderValue - lastFine) * incrFine;
  value = inputValue + incr;
  lastFine = sliderValue;
  value = max ( 0.0f, min ( maxValue, value ) );
  setFloatValue ( input, value );
}

//////////////////////
// callbacks
/////////////////////

//
// slider callbacks
//

#ifdef NEW_STRING_MODEL

static void sliderHangerMassCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 100.0, lastFine );
  float f = getFloatValue ( input );
  theString->setTension (double(f) * 9.8);
}

static void sliderHangerMassFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 100.0, lastFine );
  float f = getFloatValue ( input );
  theString->setTension ( double(f) * 9.8 );
}


static void sliderDecayTimeCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 10.0, lastFine );
  float f = getFloatValue ( input );
  theString->setDecayTime ( double(f) );
}

static void sliderDecayTimeFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 10.0, lastFine );
  float f = getFloatValue ( input );
  theString->setDecayTime ( double(f) );
}


static void sliderMassDensityCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 10.0, lastFine );
  float f = getFloatValue ( input );
  theString->setMassDensity (  double(f) );
}

static void sliderMassDensityFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 10.0, lastFine );
  float f = getFloatValue ( input );
  theString->setMassDensity ( double(f) );
}

#else

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


#endif
// #ifdef NEW_STRING_MODEL


static void sliderVibFreqCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, 1000.0, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorFreq = double(f);

}

static void sliderVibFreqFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, 1000.0, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorFreq = double(f);

}



static void sliderVibAmpCoarseCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleCoarse ( o, input, maxVibAmp, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorAmplitude = double(f);
}

static void sliderVibAmpFineCallback(Fl_Widget* o, void* input) 
{
  static float lastFine;
  handleFine ( o, input, maxVibAmp, lastFine );
  float f = getFloatValue ( input );
  theString->vibratorAmplitude = double(f);
}

//
// Fl_Input callbacks
//

#ifdef NEW_STRING_MODEL

static void inputHangerMassCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->setTension ( double(f) * 9.8 ); // tension = mg
}

static void inputDecayTimeCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->setDecayTime ( double(f) ); 
}

static void inputMassDensityCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->setMassDensity ( double(f) );
}

#else

static void inputTensionCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->Ktension = double(f);
}

static void inputDampingCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->Kdamping = double(f);
}

#endif
// #ifdef NEW_STRING_MODEL

static void inputVibFreqCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->vibratorFreq = double(f);
}

static void inputVibAmpCallback ( Fl_Widget* o, void *theFloat )
{
  float f = getFloatValue ( o );
  theString->vibratorAmplitude = double(f);
}


// 
// button callbacks
//

static void vibOnCallback ( Fl_Widget* o )
{
  theString->vibratorPower ( ((Fl_Light_Button*) o)->value() );
}

static void vibConstantPowerCallback ( Fl_Widget* o )
{
  theString->vibratorConstantPower = ((Fl_Light_Button*) o)->value();
}

static void sineButtonCallback ( Fl_Widget *o )
{
  theString->setVibratorWaveform ( 0 );
}

static void sawtoothButtonCallback ( Fl_Widget *o )
{
  theString->setVibratorWaveform ( 1 );
}


/////////////////////////////////////////////////////////////

Fl_Group *
makeCoarseFineControl ( int w, int h, const char *label, 
			Fl_Callback *fCoarseSlider, 
			Fl_Callback *fFineSlider,
			Fl_Callback *fInput,
			double initValue )
{
  int inputHeight = 20;
  int inputWidth = 50;
  int horizMargin = 30;
  //  int horizOffset = w/5;
  int horizRemain = w - inputWidth - 2*horizMargin;
  int horizRemainEnd = horizRemain - horizMargin;

  Fl_Group *group = new Fl_Group ( 0,0, w, h, NULL );


  Fl_Float_Input *input = new Fl_Float_Input ( horizMargin,0, inputWidth, inputHeight, label );
  input->labelsize(10);
  input->callback ( fInput );
  char buf[200];
  sprintf ( buf, floatFormat, initValue );
  input->value ( buf );
  group->add(input);

  Fl_Pack *sliderPack = new Fl_Pack ( inputWidth+horizMargin, 0, horizRemain,h );
  Fl_Hor_Nice_Slider *coarseSlider = new Fl_Hor_Nice_Slider( 0,0, horizRemainEnd, 2*h/3 );
  coarseSlider->color ( FL_GRAY, FL_DARK_YELLOW );
  coarseSlider->callback ( fCoarseSlider, input );

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
  o->callback ( vibConstantPowerCallback );
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
					  inputVibFreqCallback,
					  initVibFreq);
  widgetPacker->add(pack3);
  Fl_Group *pack4 = makeCoarseFineControl(width,coarsefineHeight,"Vib. Amp.", 
					  sliderVibAmpCoarseCallback, 
					  sliderVibAmpFineCallback, 
					  inputVibAmpCallback,
					  initVibAmp);
  widgetPacker->add(pack4);
  widgetPacker->end();

  return widgetPacker;
}

#ifdef NEW_STRING_MODEL

Fl_Group* 
makeNewStringControls(int x, int y, int width, int height, int coarsefineHeight )
{
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  Fl_Group *pack1 = makeCoarseFineControl ( width,coarsefineHeight,"Hanger Mass", 
					    sliderHangerMassCoarseCallback, 
					    sliderHangerMassFineCallback, 
					    inputHangerMassCallback,
					    initHangerMass);

  widgetPacker->add(pack1);

  Fl_Group *pack2 = makeCoarseFineControl ( width,coarsefineHeight,"Decay Time", 
					    sliderDecayTimeCoarseCallback, 
					    sliderDecayTimeFineCallback, 
					    inputDecayTimeCallback,
					    initDecayTime);
  widgetPacker->add(pack2);

  Fl_Group *pack3 = makeCoarseFineControl ( width,coarsefineHeight,"Mass Density", 
					    sliderMassDensityCoarseCallback, 
					    sliderMassDensityFineCallback, 
					    inputMassDensityCallback,
					    initMassDensity);
  widgetPacker->add(pack3);

  widgetPacker->end();

  return widgetPacker;
}

#else

Fl_Group* 
makeStringControls(int x, int y, int width, int height, int coarsefineHeight )
{
  //  Fl_Pack *widgetPacker = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  Fl_Pack *widgetPacker = new Fl_Pack( x,y,width,height);
  widgetPacker->spacing(10);

  Fl_Group *pack1 = makeCoarseFineControl ( width,coarsefineHeight,"Tension", 
					    sliderTensionCoarseCallback, 
					    sliderTensionFineCallback, 
					    inputTensionCallback,
					    initTension);

  widgetPacker->add(pack1);

  Fl_Group *pack2 = makeCoarseFineControl ( width,coarsefineHeight,"Damping", 
					    sliderDampingCoarseCallback, 
					    sliderDampingFineCallback, 
					    inputDampingCallback,
					    initDamping);
  widgetPacker->add(pack2);
  widgetPacker->end();

  return widgetPacker;
}
#endif
// #ifdef NEW_STRING_MODEL


Fl_Pack* 
makeApparatusControls ( int winWidth, int winHeight, int offsetWidgets, int coarsefineHeight )
{
  Fl_Pack *packer = new Fl_Pack( 0, offsetWidgets, winWidth, winHeight - offsetWidgets );
  packer->add(makeVibControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets, coarsefineHeight));

#ifdef NEW_STRING_MODEL
  packer->add(makeNewStringControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets, coarsefineHeight));
#else
  packer->add(makeStringControls(0,offsetWidgets,winWidth/2,winHeight/2-offsetWidgets, coarsefineHeight));
#endif

  packer->add ( new Fl_Window ( winWidth/2, winHeight/2 ) ); // ?????
  packer->end();
  return packer;
}
