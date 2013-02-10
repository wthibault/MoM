//
// physical string model - 
//

#pragma once
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <pthread.h>
#include <RtAudio.h>
#include <fftw3.h>

class RingBuffer {
public:
  RingBuffer ( unsigned int numFrames ) {
    totalFrames = 2*numFrames;
    data = new double [ totalFrames ];
    tail = 0;
  }
  ~RingBuffer() { delete data; }
  inline void append ( double *buffer, unsigned int nFrames ) {
    double *d = data + tail;
    if ( d+2*nFrames > data+totalFrames ) {
      // just fix tail each time
      for ( int i = 0; i < 2 * nFrames; i++ ) {
	*d++ = *buffer++;
	tail = (tail + 1) % totalFrames;
      }
    } else {
      for ( int i = 0; i < 2 * nFrames; i++ ) {
	*d++ = *buffer++;
      }
      tail = (tail + 2*nFrames) % totalFrames;
    }
  }
  inline void copyReversed ( double *dest, unsigned int numFrames, int stride ) {
    int index = tail;
    for ( int i = 0; i < numFrames; i++ ) {
      *dest++ = data[index];
      index = (index - stride);
      if (index < 0)
	index += totalFrames;
    }
  }
  int tail;
  int totalFrames;
  double *data;
};


//////////////////////////////////////////////////////


struct StringModel {
  StringModel ( int n, 
		double _Ktension, 
		double _Kdamping, 
		int _stepspersample,
		int _sampleRate,
		int _bufferFrames );
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

  inline double curvature2 ( int i, double *y );
  inline void updateElement2 ( int i );
  inline void updateElement1 ( int i );
  void computeSamples ( double *outputBuffer, unsigned int nBufferFrames );

  void analyze(double *buffer, unsigned int nBufferFrames);

  double levelDetect ( double *buffer, unsigned int nFrames );
  void compress ( double *soundout, unsigned int nBufferFrames );

  int simulationStepsPerSample;
  int numMasses;
  double Ktension;
  double Kdamping;
  int sampleRate;
  int bufferFrames;
  double *y;
  double *yold;
  double *yolder;
  double *v;
  unsigned int seed;
  pthread_mutex_t lock;
  bool   vibratorOn;
  double vibratorFreq;
  double vibratorAmplitude;
  double vibratorPhase;

  double compressionThreshold;
  double compressionRatio;

  fftw_plan     fftwPlan;
  double       *fftwIn;
  fftw_complex *fftwOut;
  int           numFramesToAnalyze;
  RingBuffer    ringBuffer;
};

