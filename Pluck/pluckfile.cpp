#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

int main ( int argc, char **argv )
{
  float freq = 440.0;  // Hz
  float dur = 1.0;     // sec
  float sr = 44100;    // sample rate
  int   i;

  for ( i = 1; i < argc; i++ ) {
    if ( strcmp ( argv[i], "-freq" ) == 0 ) {
      sscanf ( argv[i+1], "%f", &freq );
      i++;
    } else if ( strcmp ( argv[i], "-dur" ) == 0 ) {
      sscanf ( argv[i+1], "%f", &dur );
      i++;
    }
  }

  float period = 1.0 / freq;
  // ignores the 1/2 sample delay added by the lowpass
  int   delayInSamples = int(period * sr); 
  float del[delayInSamples];

  long seed = 666;
  srand48 ( seed );
  // initialize with random values;
  for ( i = 0; i < delayInSamples; i++ ) {
    del[i] = float( drand48() ) * 2.0 - 1.0;
  }

  int pos = 0;
  float ydel = 0.0;
  float recirc;
  float y;
  // run for duration
  for ( i = 0; i < dur * sr; i++ ) {
    y = del[pos];
    recirc = (y + ydel) / 2.0;
    del[pos] = recirc;
    ydel = y;
    pos++;
    if ( pos >= delayInSamples )
      pos = 0;
    short samp = short ( y * 32767.0 );
    write ( 1, &samp, 2 );
  }
}
