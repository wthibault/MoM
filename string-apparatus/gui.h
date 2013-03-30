#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Pack.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Hor_Nice_Slider.H>
#include <FL/Fl_Light_Button.H>
#include <FL/Fl_Toggle_Button.H>
#include <map>
#include <string>

Fl_Pack *makeApparatusControls ( int w, int h, int offsetWidgets, int coarsefineHeight );

void setSlider ( const char* sliderName, float value );

extern std::map<std::string, Fl_Slider*> guiSliders;
