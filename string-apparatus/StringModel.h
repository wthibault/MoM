//
// physical string model - 
//

#pragma once
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <pthread.h>
#include <RtAudio.h>

struct StringModel {
  StringModel ( int n, 
		double _Ktension, 
		double _Kdamping, 
		int _stepspersample );
  ~StringModel();

  void print();
  void reset();
  void pluck();
  void pluckvel();
  void toggleVibrator();
  inline double linearToDecibels ( double amp );
  inline double decibelsToLinear ( double db );
  inline void clip ( double *s );

  // for RtAudio
  static int audioCallback ( void *outputBuffer, 
			     void *inputBuffer, 
			     unsigned int nBufferFrames,
			     double streamTime, 
			     RtAudioStreamStatus status, 
			     void *userData );

  void computeSamples ( double *outputBuffer, unsigned int nBufferFrames );

  double levelDetect ( double *buffer, unsigned int nFrames );
  void compress ( double *soundout, unsigned int nBufferFrames );

  int simulationStepsPerSample;
  int numMasses;
  double Ktension;
  double Kdamping;
  int sampleRate;
  double *y, *yold, *v;
  unsigned int seed;
  pthread_mutex_t lock;
  bool   vibratorOn;
  double vibratorFreq;
  double vibratorAmplitude;
  double vibratorPhase;

  double compressionThreshold;
  double compressionRatio;

};
