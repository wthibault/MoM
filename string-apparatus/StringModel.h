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

////////////////////////////////////////////////////



class Histogram {
 public:
  Histogram ( unsigned int num_bins, double minRange, double maxRange ) 
    : numBins ( num_bins ),
    minVal ( minRange ),
    maxVal ( maxRange ),
    bins ( new unsigned int[num_bins] ) 
      {
	bins = new unsigned int [numBins];
	clear();
      }
  Histogram () 
    :  numBins ( 256 ), 
       minVal ( -1.0 ), 
       maxVal ( 1.0 ) 
	 { 
	   bins = new unsigned int [numBins];
	   clear(); 
	 }
  ~Histogram () { delete bins; }
  void clear () {
    for (int i = 0; i < numBins; i++ ) {
      bins[i] = 0;
    }
  }
  inline void update ( double val ) {
    if ( val < minVal ) minVal = val;
    if ( val > maxVal ) maxVal = val;
#if USE_HISTOGRAM_BINS
    int bin;
    if ( maxVal - minVal < 1e-6 )
      bin = 0;
    else
      bin = static_cast<int> ( (numBins-1) * ( (val - minVal) / (maxVal - minVal) ) ); // truncate
    bins[bin]++;
#endif
  }
  unsigned int *bins;
  unsigned int numBins;
  double minVal;
  double maxVal;
};



////////////////////////////////////////////////////



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
	data[tail] = *buffer++;
	tail = (tail + 1) % totalFrames;
      }
    } else {
      for ( int i = 0; i < 2 * nFrames; i++ ) {
	*d++ = *buffer++;
      }
      tail = (tail + 2*nFrames); // % totalFrames;
    }
  }

  inline double windowFunction ( double x ) {
    // x should be between 0 and 1
    // type this into octave to visualize:
    // function y=g(x)
    // y= 0.5 .* ( 1 .+ cos ( pi .* (x*2-1) ) );
    // endfunction 
    // plot(x,g(x))

    // Hanning window, definite integral over [0,1] is 1.0
    return 0.5 * ( 1.0 + cos ( M_PI * (x*2-1) ) ); 
  }

  inline void copyReversed ( double *dest, unsigned int numFrames, int stride, bool windowed = true ) {
    int index = tail;
    if ( windowed ) {
      double incr = 1.0 / numFrames;
      for ( int i = 0; i < numFrames; i++ ) {
	//	*dest++ = windowFunction ( index * incr ) * data[index];
	*dest++ = windowFunction ( i * incr ) * data[index];
	index = (index - stride);
	if (index < 0)
	  index += totalFrames;
      }
    } else {
      for ( int i = 0; i < numFrames; i++ ) {
	*dest++ = data[index];
	index = (index - stride);
	if (index < 0)
	  index += totalFrames;
      }
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

  Histogram    *histograms;
};

