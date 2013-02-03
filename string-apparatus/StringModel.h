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
  StringModel ( int n, float _Ktension, float _Kdamping, int _stepspersample )
    : numMasses(n), 
      Ktension(_Ktension),
      Kdamping(_Kdamping),
    sampleRate ( 44100 ),
      simulationStepsPerSample(_stepspersample),
      vibratorOn ( false ),
      vibratorFreq ( 100.0f ),
      vibratorAmplitude ( 0.001f ),
      vibratorPhase ( 0.0 ),
      compressionThreshold ( -10.0 ),
      compressionRatio ( 0.5 )
  {
    std::cout << "StringModel N,K_t,K_d,ss: " 
	      << numMasses << ',' 
	      << Ktension  << ',' 
	      << Kdamping  << ',' 
	      << simulationStepsPerSample
	      << std::endl;

    // allocate displacement and velocity arrays
    y = new float[numMasses];
    yold = new float[numMasses];
    v = new float[numMasses];
    // initialize displacements and velocities
    for (int i = 0; i < numMasses; i++ ) {
      v[i]  = 0.0f;
      yold[i] = y[i] = 0.0f;
    }
    // init the mutex
    pthread_mutex_init ( &lock, NULL );
    // init the RNG
    seed = (unsigned int) time(NULL);
    // pluck it
    pluck();
  }

  ~StringModel() { delete y; delete yold, delete v; }

  void print() {
    pthread_mutex_lock ( &lock );
    std::cout << "N,K_t,K_d,ss: " 
	      << numMasses << Ktension << Kdamping << simulationStepsPerSample
	      << std::endl;
    for (int i = 0; i < numMasses; i++) {
      std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    pthread_mutex_unlock ( &lock );
  }


  void reset() {
    pthread_mutex_lock ( &lock );
    for (int i = 0 ; i < numMasses; i++ ) {
      v[i] = y[i] = yold[i] = 0.0;
    }
    pthread_mutex_unlock ( &lock );
  }


  void pluck() {
    int pluckAt = float(rand_r(&seed) / float(RAND_MAX)) * (numMasses-2) + 1;
    //    std::cout << pluckAt << std::endl;
    pthread_mutex_lock( &lock );
    float maxDisp = 0.001;
    float upSlope = maxDisp / pluckAt;
    float downSlope = maxDisp / (numMasses-pluckAt);
    for (int i = 1; i< numMasses-2;i++ ) {
      if (i <= pluckAt )
	yold[i] = i*upSlope;
      else
	yold[i] = maxDisp - (i-pluckAt)*downSlope; 
    }
    pthread_mutex_unlock ( &lock );
  }

  void pluckvel() {
    int pluckAt = float(rand_r(&seed) / float(RAND_MAX)) * (numMasses-2) + 1;
    //    std::cout << pluckAt << std::endl;
    pthread_mutex_lock( &lock );
    v[pluckAt] = 3e-4 * Ktension;
    pthread_mutex_unlock ( &lock );
  }

  void toggleVibrator() {
    vibratorOn = !vibratorOn;
  }

  inline float linearToDecibels ( float amp ) {
    return 20.0f * log10 ( amp );
  }

  inline float decibelsToLinear ( float db ) {
    return pow ( 10.0, (0.05 * db) );
  }

  inline void clip ( float *s ) {
    if (fabs(*s) > 1.0) {
      *s = *s < 0.0 ? -1.0 : 1.0;
    }
  }

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
