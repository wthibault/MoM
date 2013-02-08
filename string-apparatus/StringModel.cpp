// 
// StringModel.cpp
//

#include "StringModel.h"

// BEGIN pure c++ audio computation
// --no calls to the outside world allowed



// simple lerper class
class Parameter {
  Parameter() : targetValue(0.0f),
		previousValue ( 0.0f )
  {}
  inline double getValueAt(double u) {
    return lastValue = u * targetValue + (1.0f-u) * previousValue;
  }
  inline void setTargetValue( double target ) {
    targetValue = target;
    previousValue = lastValue;
  }
  inline void setCurrentValue ( double v ) {
    previousValue = lastValue = v;
  }
  double targetValue;
  double previousValue;
  double lastValue;
};


///
/// computeSamples creates a buffer's worth of output samples
///
void
StringModel::computeSamples ( double *soundout, unsigned int nBufferFrames )
{


  int iters = nBufferFrames * simulationStepsPerSample;
  double *buf = soundout;

  for ( int t = 0; t < iters; t++ ) {

    double sum = 0;
    int n = numMasses;
    int i;
    double accel;
    if ( vibratorOn ) {
      // XXX interpolate params
      y[0] = vibratorAmplitude * sin ( vibratorPhase );
      vibratorPhase += (vibratorFreq / sampleRate) / simulationStepsPerSample;
      while ( vibratorPhase > 2*M_PI )
	vibratorPhase -= 2*M_PI;
    } else {
      y[0] = 0.0f;
    }

    // XXX omp not working : #pragma omp parallel for reduction(+:sum) private(accel)
    for ( i = 1; i < n-1; i++ ) {
	accel = Ktension * (yold[i+1] + yold[i-1] - 2*yold[i]);
	v[i] += accel;
	v[i] *= Kdamping;
	y[i] = yold[i] + v[i];
	sum = sum + y[i];
    }

    double *tmp = y;
    y = yold;
    yold = tmp;

    if ( t % simulationStepsPerSample == 0 ) {


#if 0
      // summed output
      *buf++ = sum;
      *buf++ = sum;
#else
      // per channel pickup placement
      *buf++ = decibelsToLinear(30.0) * yold[numMasses/13];
      *buf++ = decibelsToLinear(30.0) * yold[numMasses/5];
#endif
    }

  }

  //
  // Limit/compress (dynamic volume adjustment)
  //

  compress ( soundout, nBufferFrames );

}

void
StringModel::compress ( double *soundout, unsigned int nBufferFrames )
{
  // exponential smoothing for envelope follower
  static double rmsLevel = 0.0;
  double rmsLevelMeasurement = levelDetect ( soundout, nBufferFrames );
  double alpha = 0.1;// XXX
  double beta = 0.1;// XXX
  if ( rmsLevelMeasurement < rmsLevel )
    rmsLevel = beta * rmsLevelMeasurement + (1-beta) * rmsLevel;
  else
    rmsLevel = alpha * rmsLevelMeasurement + (1-alpha) * rmsLevel;
  int numFrames = nBufferFrames;
  double *sample = soundout;
  double gain;

  if ( rmsLevel > compressionThreshold ) {
    gain = decibelsToLinear ( - (rmsLevel - compressionThreshold) * compressionRatio );
  } else {
    gain = 1.0f;
  }

  while (numFrames--) {
    clip(sample);
    *sample++ = *sample * gain;
    clip(sample);
    *sample++ = *sample * gain;
  }
}

inline
double
StringModel::levelDetect ( double *buffer, unsigned int nFrames )
{
  // compute RMS average over the buffer
  double sumOfSquares = 0.0f;
  double val;
  for (unsigned int i = 0; i < nFrames*2; i++) {
    val = *buffer++;
    sumOfSquares += val * val;
  }
  return linearToDecibels (sqrt ( sumOfSquares / nFrames*2 ) );
}







//
// END pure c++ callback section
//

//////////////////////////////////////////////////////////////////////////////////////////

//
// the RtAudio callback
//

int
StringModel::audioCallback ( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	 double streamTime, RtAudioStreamStatus status, void *userData )
{

  if (status == RTAUDIO_INPUT_OVERFLOW)
    std::cout << "string -- overflow" << std::endl;

  double *soundout = (double *)outputBuffer;
  StringModel *s = (StringModel *)userData;

  pthread_mutex_lock(&(s->lock));

  s->computeSamples ( static_cast<double*>(outputBuffer), nBufferFrames );
  
  pthread_mutex_unlock ( &(s->lock) );
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

StringModel::StringModel ( int n, double _Ktension, double _Kdamping, int _stepspersample )
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
  y = new double[numMasses];
  yold = new double[numMasses];
  v = new double[numMasses];
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

StringModel::~StringModel() 
{ 
  delete y; 
  delete yold; 
  delete v; 
}

void 
StringModel::print() 
{
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


void 
StringModel::reset() 
{
  pthread_mutex_lock ( &lock );
  for (int i = 0 ; i < numMasses; i++ ) {
    v[i] = y[i] = yold[i] = 0.0;
  }
  pthread_mutex_unlock ( &lock );
}


void 
StringModel::pluck() 
{
  int pluckAt = double(rand_r(&seed) / double(RAND_MAX)) * (numMasses-2) + 1;
  //    std::cout << pluckAt << std::endl;
  pthread_mutex_lock( &lock );
  double maxDisp = 0.001;
  double upSlope = maxDisp / pluckAt;
  double downSlope = maxDisp / (numMasses-pluckAt);
  for (int i = 1; i< numMasses-2;i++ ) {
    if (i <= pluckAt )
      yold[i] = i*upSlope;
    else
      yold[i] = maxDisp - (i-pluckAt)*downSlope; 
  }
  pthread_mutex_unlock ( &lock );
}

void 
StringModel::pluckvel() 
{
  int pluckAt = double(rand_r(&seed) / double(RAND_MAX)) * (numMasses-2) + 1;
  //    std::cout << pluckAt << std::endl;
  pthread_mutex_lock( &lock );
  v[pluckAt] = 3e-4 * Ktension;
  pthread_mutex_unlock ( &lock );
}

void 
StringModel::toggleVibrator() {
  vibratorOn = !vibratorOn;
}

inline double 
StringModel::linearToDecibels ( double amp ) 
{
  return 20.0f * log10 ( amp );
}

inline double 
StringModel::decibelsToLinear ( double db ) 
{
  return pow ( 10.0, (0.05 * db) );
}

inline void 
StringModel::clip ( double *s ) 
{
  if (fabs(*s) > 1.0) {
    *s = *s < 0.0 ? -1.0 : 1.0;
  }
}


