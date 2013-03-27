// 
// StringModel.cpp
//

#include "StringModel.h"


// XXX use a vibration modulation transfer function cf. http://spie.org/x34359.xml

//////////////////////////////////////////
// BEGIN pure c++ audio computation
// --no calls to the outside world allowed

inline
void
swap3 ( double * &a, double * &b, double * &c )
{
  double *tmp = c;
  c = b;
  b = a;
  a = tmp;
}
  


inline
void
StringModel::updateElement1 ( int i ) 
{
  // !!! dont call this for end elements!
  y[i] = c1 * ( 2*yold[i] - c2*yolder[i] + c3*(yold[i+1]+yold[i-1]-2*yold[i]) );
  //  sum = sum + y[i]; // don't need this if using pickups
}



inline
double
StringModel::adjustedVibratorAmplitude()
{
  if ( vibratorConstantPower )
    return vibratorAmplitude / (vibratorFreq*vibratorFreq);
  else
    return vibratorAmplitude;
}

///
/// computeSamples creates a buffer's worth of output samples
///
void
StringModel::computeSamples ( double *soundout, unsigned int nBufferFrames )
{


  int iters = nBufferFrames * simulationStepsPerSample;
  double *buf = soundout;

  lock();

  for ( int t = 0; t < iters; t++ ) {

    double sum = 0;
    int n = numMasses;
    int i;
    double accel;

    // XXX REFACTOR
    // see a plot:
    // octave:11> n=0:44100;
    // octave:12> function y=digitalOsc(f,sr,n)
    // > y = sin ( 2 .* pi .* f .* ( n ./ sr ));
    // > endfunction
    // octave:14> plot(n,digitalOsc(1,44100,n))

    if ( vibratorOn ) {
      if ( vibratorWaveform == 0 ) 
	// sine
	y[0] = adjustedVibratorAmplitude() * sin ( vibratorPhase );
      else
	// sawtooth
	y[0] = adjustedVibratorAmplitude() * (vibratorPhase / (2*M_PI));

      vibratorPhase += ( 2*M_PI * vibratorFreq / sampleRate) 
	                / simulationStepsPerSample;
      while ( vibratorPhase > 2*M_PI )
	vibratorPhase -= 2*M_PI;
    } else {
      y[0] = 0.0;
    }


    // update positions
    for ( i = 1; i < n-1; i++ ) {
      updateElement1 ( i );
      histograms[i].update ( y[i] );
    }

    // triple buffer needed for position update
    swap3 ( y, yold, yolder );

    if ( t % simulationStepsPerSample == 0 ) {


#if 0
      // summed output
      *buf++ = sum;
      *buf++ = sum;
#else
      // per channel pickup placement
      *buf++ = decibelsToLinear(30.0) * yold[numMasses/13]; // XXX setPickupPositions (l,r)
      *buf++ = decibelsToLinear(30.0) * yold[numMasses/5];
#endif
    }

  }


  freshHistograms = true;

  unlock();


  // 
  // ship a copy off for FFT analysis
  //

  analyze ( soundout, nBufferFrames );

  //
  // Limit/compress (dynamic volume adjustment)
  //

  compress ( soundout, nBufferFrames );

}

void
StringModel::analyze (double *buffer, unsigned int nBufferFrames)
{
  // dont do the analysis here, but tuck away this buffer for later
  // analysis by FFTPrimitive, as needed for rendering.
  ringBuffer.append ( buffer, nBufferFrames );
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

  s->computeSamples ( static_cast<double*>(outputBuffer), nBufferFrames );

  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


StringModel::StringModel ( int n, 
			   double _Ktension, 
			   double _massDensity,
			   double _decayTime,
			   double _length,
			   int _stepspersample,
			   int _sampleRate,
			   int _bufferFrames )
    : numMasses(n), 
      Ktension(_Ktension),
      massDensity ( _massDensity ),
      decayTime ( _decayTime ),
      length ( _length ),
      sampleRate ( _sampleRate ),
      bufferFrames ( _bufferFrames ),
      simulationStepsPerSample(_stepspersample),
      vibratorOn ( false ),
      vibratorWaveform ( 0 ),
      vibratorFreq ( 100.0f ),
      vibratorAmplitude ( 0.001f ),
      vibratorPhase ( 0.0 ),
      vibratorConstantPower ( false ),
      compressionThreshold ( -10.0 ),
      compressionRatio ( 0.75 ),
      numFramesToAnalyze ( 8 * _bufferFrames ),
      ringBuffer ( numFramesToAnalyze * 2 ),
      freshHistograms ( false )
{
  std::cout << "StringModel N,K_t,mu,tau,ss: " 
	    << numMasses << ',' 
	    << Ktension  << ',' 
	    << massDensity  << ',' 
	    << decayTime  << ',' 
	    << simulationStepsPerSample
	    << std::endl;

  // allocate displacement and velocity arrays
  y = new double[numMasses];
  yold = new double[numMasses];
  yolder = new double[numMasses];

  // initialize displacements and velocities
  for (int i = 0; i < numMasses; i++ ) {
    yolder[i] = yold[i] = y[i] = 0.0f;
  }

  // precompute constants used in update
  precomputeConstants();
  
  // init the mutex
  initLock();

  // init the RNG
#ifdef __WIN32__
#else  
  seed = (unsigned int) time(NULL);
#endif

  // pluck it
  //  pluck();

  // set up the fft
  // XXX move this into FFTPrimitive XXX ???
  fftwIn = (double *)        fftw_malloc(sizeof(double)*numFramesToAnalyze);
  fftwOut = (fftw_complex *) fftw_malloc ( sizeof(fftw_complex) * (numFramesToAnalyze/2 + 1));
  fftwPlan = fftw_plan_dft_r2c_1d ( numFramesToAnalyze, fftwIn, fftwOut, FFTW_ESTIMATE );

  // setup the histograms: one per mass! for "vibration modulation transfer function" estimate
  //  for ( int i = 0; i < numMasses; i++ ) {
    histograms = new Histogram [numMasses]; // ( 256, -1.0, 1.0 );
    //  }

}

const double EPSILON = 1e-4;

void
StringModel::precomputeConstants() {
  double dt = 1.0 / (sampleRate * simulationStepsPerSample);
  double dx = length / numMasses;
  c1 = 1.0 / (1.0 + (dt/(2*decayTime)));
  c2 = 1.0 - dt/(2*decayTime);
  double c = sqrt(Ktension / massDensity);
  c3 = c * dt / dx;
  c3 *= c3;
  c3 = std::max ( EPSILON, std::min ( 1.0, c3 ) );
  std::cout << "c1,c2,c3 = " << c1 << ' ' << c2 << ' ' << c3 << std::endl;
}


double 
clamp(double val)
{
  if ( val <= EPSILON) return EPSILON;
  else return val;
}

void
StringModel::setTension ( double T )
{
  T = clamp(T);
  Ktension = T;
  precomputeConstants();
}

void
StringModel::setMassDensity ( double mu )
{
  //  std::cout << "setMassDensity mu: " << mu << std::endl;
  //  const double EPSILON = 1e-2;
  //  if ( std::abs(mu) <= EPSILON) mu = EPSILON;
  massDensity = clamp(mu);
  precomputeConstants();
}

void
StringModel::setDecayTime ( double tau )
{
  //  if ( tau <= 0 ) return;
  decayTime = clamp(tau);
  precomputeConstants();
}

void
StringModel::setVibratorWaveform ( int wave )
{
  vibratorWaveform = wave;
}



StringModel::~StringModel() 
{ 
  delete y; 
  delete yold; 
  delete yolder; 

  fftw_destroy_plan(fftwPlan);
  fftw_free(fftwIn );
  fftw_free(fftwOut);
}

void 
StringModel::print() 
{
  lock();
  std::cout << "N,K_t,mu,tau,ss: " 
	    << numMasses << Ktension << massDensity << decayTime << simulationStepsPerSample
	    << std::endl;
  for (int i = 0; i < numMasses; i++) {
    std::cout << y[i] << " ";
  }
  std::cout << std::endl;
  unlock();
}

void 
StringModel::reset() 
{
  lock();
  for (int i = 0 ; i < numMasses; i++ ) {
    y[i] = yold[i] = yolder[i] = 0.0;
  }
  unlock();
}


void 
StringModel::pluck() 
{
  int pluckAt;
#  ifdef __WIN32__
    pluckAt = double(rand() / double(RAND_MAX)) * (numMasses-2) + 1;
#else
    pluckAt = double(rand_r(&seed) / double(RAND_MAX)) * (numMasses-2) + 1;
#endif
  //    std::cout << pluckAt << std::endl;
    lock();
  double maxDisp = 0.001;
  double upSlope = maxDisp / pluckAt;
  double downSlope = maxDisp / (numMasses-pluckAt);
  for (int i = 1; i< numMasses-2;i++ ) {
    if (i <= pluckAt )
      yold[i] = i*upSlope;
    else
      yold[i] = maxDisp - (i-pluckAt)*downSlope; 
  }
  unlock();
}

void 
StringModel::pluckvel() 
{
# ifdef __WIN32__
  int pluckAt = double(rand() / double(RAND_MAX)) * (numMasses-2) + 1;
# else
  int pluckAt = double(rand_r(&seed) / double(RAND_MAX)) * (numMasses-2) + 1;
# endif
  //    std::cout << pluckAt << std::endl;
  lock();
  // XXX NOOP
  unlock();
}

void 
StringModel::toggleVibrator() {
  vibratorOn = !vibratorOn;
}

void 
StringModel::vibratorPower(bool on) {
  vibratorOn = on;
}

inline double 
StringModel::linearToDecibels ( double amp ) 
{
  return 20.0 * log10 ( amp );
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


