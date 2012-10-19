#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


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

  // string displacements and velocities
  float *y, *yold;
  float *v, *vold;

  //  const int NUM_MASSES = 1000;
  //  const float Ktension = 0.1;
  //  const float Kdamping = 0.9;
  //  const float duration = 10.0;

  if ( argc != 5 ) {
    std::cout << "usage: " << argv[0] << " NUM_MASSES durationInSecs Ktension Kdamping\n";
    return 2;
  }
  int NUM_MASSES  = atoi ( argv[1] );
  float duration = atof ( argv[2] );
  float Ktension = atof ( argv[3] );
  float Kdamping = atof ( argv[4] );
  const int PICKUP_POS = NUM_MASSES / 7;    // wtf

  // open output file
  FILE *f = fopen ( "string.raw", "wb" );
  if (!f) {
    std::cout << "can't open output file\n";
    return 1;
  }
  // allocate displacement and velocity arrays
  y = new float[NUM_MASSES];
  yold = new float[NUM_MASSES];
  v = new float[NUM_MASSES];

  // initialize displacements (pluck it!) and velocities
  for (int i = 0; i < NUM_MASSES; i++ ) {
    v[i]  = 0.0f;
    yold[i] = y[i] = 0.0f;
    if (i == NUM_MASSES/2 )
      yold[i] = 1.0; // impulse at string center
  }

  // run simulation for desired period:
  int numIters = duration * 44100; 
  //#pragma omp parallel for
  for ( int t = 0; t < numIters; t++ ) {
    // for each mass element
    float sum = 0;
    for ( int i = 0; i < NUM_MASSES; i++ ) {
      //   if boundary element
      //      handle boundary element
      if ( i == 0 || i == NUM_MASSES-1 ) {
      } else {
	//   else
	//      compute acceleration as scaled sum of differences with neighbors
	float accel = Ktension * (yold[i+1] + yold[i-1] - 2*yold[i]);
	//      add accel to velocity
	v[i] += accel;
	v[i] *= Kdamping;
	//      add velocity to position
	y[i] = yold[i] + v[i];
	sum += y[i];
      }
    }
    //   swap displacement buffers
    float *tmp = y;
    y = yold;
    yold = tmp;

    //   output a sample as sum of displacements
    fwrite ( &sum, sizeof(float), 1, f );

    // output a sample at pickup
    //fwrite ( &y[PICKUP_POS], sizeof(float), 1, f );
  }
  // close output file
  fclose ( f );
}
