//
// ui.cpp
//

#include "ui.h"
#include <FTGL/ftgl.h>
#include <iostream>
#include <cmath>
#include "StringModel.h"






  ////////////////////////////////////////////////////////////////////





UI::UI ( const char *path )
  : font ( path )
{
  if (font.Error()) {
    std::cout << "problem loading font " << path << std::endl;
  }  
}

UI::~UI()
{
}

void
UI::printAt ( float x, float y, std::string msg )
{
  font.Render(msg.c_str(), -1, FTPoint(x,y));
}


void
UI::draw (int width, int height, StringModel *str)
{
  float fProportion = 20;
  float fSize = fmin(width,height) / fProportion;
  font.FaceSize(fSize);
  char buf[100];

  sprintf(buf,"tension   %9.7lf", str->Ktension);
  printAt ( 0, 0.9*height, buf );

  sprintf(buf,"damping  %9.7lf", str->Kdamping);
  printAt ( 0, 0.8*height, buf );

  sprintf(buf,"vib freq    %9.4lf Hz", str->vibratorFreq);
  printAt ( 0, 0.7*height, buf );
  
  sprintf(buf,"vib amp    %9.4lf mm", str->vibratorAmplitude*1000);
  printAt ( 0, 0.6*height, buf );

  if ( str->vibratorOn ) 
    sprintf(buf,"vib on");
  else
    sprintf(buf,"vib off");
  printAt ( 0, 0.5*height, buf );

  
}

