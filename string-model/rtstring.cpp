#include <iostream>
#include <vector>
#include <stdio.h>
#include <string>
#include <cstdlib>
#include <pthread.h>
#include <omp.h>
#include <cmath>

#include "RtAudio.h"

struct StringState {
  StringState ( int n, float _Ktension, float _Kdamping, int _stepspersample )
    : numMasses(n), 
      Ktension(_Ktension),
      Kdamping(_Kdamping),
      simulationStepsPerSample(_stepspersample)
      
  {
    std::cout << "N,K_t,K_d,ss: " 
	      << numMasses << Ktension << Kdamping << simulationStepsPerSample 
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
  ~StringState() { delete y; delete yold, delete v; }

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
      //      if (i == numMasses/2 )
      //	yold[i] = 0.5; // impulse at string center
    //  impluse at string center
    //    yold[numMasses/2] = 0.25;
    int pluckAt = float(rand_r(&seed) / float(RAND_MAX)) * (numMasses-2) + 1;
    std::cout << pluckAt << std::endl;
    pthread_mutex_lock( &lock );
    yold[pluckAt] = 0.1;
    //    v[pluckAt] = 5e-3;
    pthread_mutex_unlock ( &lock );
  }

  int simulationStepsPerSample;
  int numMasses;
  float Ktension;
  float Kdamping;
  float *y, *yold, *v;
  unsigned int seed;
  pthread_mutex_t lock;
};


int
string ( void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
	 double streamTime, RtAudioStreamStatus status, void *userData )
{

  if (status == RTAUDIO_INPUT_OVERFLOW)
    std::cout << "string -- overflow" << std::endl;

  float *soundout = (float *)outputBuffer;
  StringState *s = (StringState *)userData;

  pthread_mutex_lock(&(s->lock));




#if 0
  int nthreads, tid;
#pragma omp parallel private(nthreads, tid)
  {

  /* Obtain thread number */
  tid = omp_get_thread_num();
  printf("Hello World from thread = %d\n", tid);

  /* Only master thread does this */
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }

  }  /* All threads join master thread and disband */
#endif





  int iters = nBufferFrames * s->simulationStepsPerSample;
  for ( int t = 0; t < iters; t++ ) {
    // for each mass element
    float sum = 0;
    int n = s->numMasses;
    int i;
    float accel;
    //#pragma omp parallel for reduction(+:sum) private(accel)
    for ( i = 0; i < n; i++ ) {
      //   if boundary element
      //      handle boundary element
      if ( i == 0 || i == s->numMasses-1 ) {
      } else {
	//   else
	//      compute acceleration as scaled sum of differences with neighbors
	accel = s->Ktension * (s->yold[i+1] + s->yold[i-1] 
				     - 2*s->yold[i]);
	//      add accel to velocity
	s->v[i] += accel;
	s->v[i] *= s->Kdamping;
	//      add velocity to position
	s->y[i] = s->yold[i] + s->v[i];
	sum = sum + s->y[i];
      }
    }
    //   swap displacement buffers
    float *tmp = s->y;
    s->y = s->yold;
    s->yold = tmp;

    if ( t % s->simulationStepsPerSample == 0 ) {
      if (fabs(sum) > 1.0) {
	sum = 1.0;
	std::cout << "! " << sum << std::endl;
      }
#if 1
      *soundout++ = sum;
      *soundout++ = sum;
#else
      *soundout++ = 10.0f * s->yold[s->numMasses/13];
      *soundout++ = 10.0f * s->yold[s->numMasses/5];
#endif
    }

  }
  
  pthread_mutex_unlock ( &(s->lock) );
  return 0;
}


int main ( int argc, char **argv )
{
  // algorithm:
  // open output file
  // allocate displacement arrays and velocities
  // initialize displacements (pluck it!) and velocities
  // run simulation for desired period:
  // for each mass element
  //   if boundary element
  //      handle boundary element
  //   else
  //      compute acceleration as scaled sum of differences with neighbors
  //      add accel to velocity
  //      add velocity to position
  //   swap displacement buffers
  //   output a sample (sum of displacements)
  // close output file

  //  StringState theString ( 600, 0.75, 0.99999, 4 );
  StringState theString ( 600, 0.5, 0.99999, 4 );

  // string displacements and velocities
  float *y, *yold;
  float *v, *vold;

  // *** test the rtaudio callback
  RtAudio dac;
  if ( dac.getDeviceCount() < 1 ) {
    std::cout << "\nNo audio devices found!\n";
    exit( 0 );
  }
    
  RtAudio::StreamParameters parameters;
  parameters.deviceId = dac.getDefaultOutputDevice();
  parameters.nChannels = 2;
  parameters.firstChannel = 0;
  unsigned int sampleRate = 44100;
  unsigned int bufferFrames = 256; // 256 sample frames

  try { 
    dac.openStream ( &parameters, NULL, RTAUDIO_FLOAT32,
		     sampleRate, &bufferFrames, 
		     &string, (void *)&theString );
    dac.startStream();
  }
  catch ( RtError& e ) {
    e.printMessage();
    exit(0);
  }



#if 1
  int nthreads, tid;
#pragma omp parallel private(nthreads, tid)
  {

  /* Obtain thread number */
  tid = omp_get_thread_num();
  printf("Hello World from thread = %d\n", tid);

  /* Only master thread does this */
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }

  }  /* All threads join master thread and disband */
#endif







    char input;
    std::cout << "\nPlaying ... press 'q' to quit.\n";
    std::cout << "press 'p' to pluck\n";
    std::cout << "press 't' to tighten\n";
    std::cout << "press 'l' to loosen\n";
    bool done = false;
    do {
      std::cin.get( input );
      switch ( input ) {
      case 't' : theString.Ktension *= 1.05946; break;
      case 'l' : theString.Ktension *= 0.943876; break;
      case 'p' : theString.pluck(); break;
      case 'r' : theString.reset(); break;
      case 'd' : theString.print(); break;
      case 'q' : done = true; break;
      }
    } while (!done);
    
    try {
        // Stop the stream
        dac.stopStream();
    }
    catch (RtError& e) {
        e.printMessage();
    }

  if (dac.isStreamOpen()) dac.closeStream();

  return 0;
}
