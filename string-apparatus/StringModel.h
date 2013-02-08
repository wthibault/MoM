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
  StringModel ( int n, float _Ktension, float _Kdamping, int _stepspersample );
  ~StringModel();

  void print();
  void reset();
  void pluck();
  void pluckvel();
  void toggleVibrator();
  inline float linearToDecibels ( float amp );
  inline float decibelsToLinear ( float db );
  inline void clip ( float *s );

  static int audioCallback ( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
			   double streamTime, RtAudioStreamStatus status, void *userData );

  void computeSamples ( float *outputBuffer, unsigned int nBufferFrames );

  float levelDetect ( float *buffer, unsigned int nFrames );
  void compress ( float *soundout, unsigned int nBufferFrames );

  int simulationStepsPerSample;
  int numMasses;
  float Ktension;
  float Kdamping;
  int sampleRate;
  float *y, *yold, *v;
  unsigned int seed;
  pthread_mutex_t lock;
  bool   vibratorOn;
  float vibratorFreq;
  float vibratorAmplitude;
  float vibratorPhase;
  //  float t;

  float compressionThreshold;
  float compressionRatio;

};
