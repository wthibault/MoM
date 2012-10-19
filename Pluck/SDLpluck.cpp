#include <SDL/SDL.h>
#include <vector>
using namespace std;

class PluckedString {
public:
  int delayInSamples; // length of delay line
  float sr;    // sampleRate
  float *del;  // delay line
  int   pos;   // current pos in delay line
  float ydel;  // last output value
  PluckedString ( float freq, float sr ) {
    delayInSamples = int ( sr / freq );
    del = new float[delayInSamples];
    for ( int i = 0; i < delayInSamples; i++ ) {
      del[i] = 0.0;
    }
    ydel = 0.0;
    pos = 0;
  };
  virtual ~PluckedString () {
    delete del;
  };

  void Pluck ( float strength ) {
    const float maxAmp = 0.25;
    SDL_LockAudio();
    if ( strength > maxAmp )
      strength = maxAmp;
    for ( int i = 0; i < delayInSamples; i++ ) {
      del[i] += float ( strength * (drand48() * 2.0 - 1.0) );
      if ( del[i] < -1.0 ) del[i] = -1.0;
      if ( del[i] >  1.0 ) del[i] = 1.0;
    }
    SDL_UnlockAudio();
  };

  // gets called from the callback
  void genSamples ( short *samples, int numSamples ) {
    float y, recirc;
    for ( int i = 0; i < numSamples; i++ ) {
      y = del[pos];
      recirc = (y + ydel) / 2.0;
      del[pos] = recirc;
      ydel = y;
      pos++;
      if ( pos >= delayInSamples )
	pos = 0;

      // mix into output w/clipping
      int samp = short ( y * 32767.0 ) + *samples;
      if ( samp > 32767 ) samp = 32767;
      if ( samp < -32767 ) samp = -32767;
      *samples++ = short(samp);
    }
  };

};


vector < PluckedString* > strings;

void
initPlucker ( SDL_AudioSpec *spec )
{
  if ( spec->format != AUDIO_S16SYS ) {
    fprintf ( stderr, "unsupported format!\n");
    exit ( 3 );
  }

  float base_freq = 440.0;
  vector<float> ratios;
  ratios.push_back ( 16.0/15.0 );
  ratios.push_back ( 2.0 );
  ratios.push_back ( 5.0/3.0 );
  ratios.push_back ( 9.0/8.0 );
  ratios.push_back ( 6.0/5.0 );
  ratios.push_back ( 15.0/8.0 );
  ratios.push_back ( 16.0/7.0 );
  ratios.push_back ( 1.0 );
  ratios.push_back ( 4.0/3.0 );
  ratios.push_back ( 45.0/32.0 );
  ratios.push_back ( 3.0/2.0 );
  ratios.push_back ( 5.0/4.0 );

  for ( int i = 0; i < ratios.size(); i++ ) {
    strings.push_back (new PluckedString ( base_freq * ratios[i], 
					   spec->freq ));
  }

}

// SDL audio callback
void
plucker ( void *userdata, Uint8 *stream, int len )
{
  for ( unsigned int i = 0; i < strings.size(); i++ ) {
    strings[i]->genSamples ( (short *)stream, len/2 );
  }
}

int main ( int argc, char *argv[] ) 
{

  if ( SDL_Init ( SDL_INIT_VIDEO | SDL_INIT_TIMER | SDL_INIT_AUDIO ) < 0 ) {
    fprintf ( stderr, "cant init SDL %s\n", SDL_GetError() );
    exit ( 1 );
  }

  SDL_AudioSpec desired, obtained;
  desired.freq = 44100;
  desired.format = AUDIO_S16SYS;
  desired.channels = 1;
  desired.samples = 2048;
  desired.callback = plucker;
  desired.userdata = NULL;

  if ( SDL_OpenAudio (&desired, &obtained) < 0 ) {
    fprintf ( stderr, "Couldn't open audio: %s\n", SDL_GetError() );
    exit ( 2 );
  }

  initPlucker ( &obtained );

  SDL_Surface *screen = SDL_SetVideoMode ( 200,200,16,SDL_SWSURFACE);
  SDL_WM_SetCaption("SDLpluck", NULL);

  SDL_PauseAudio ( 0 );  // start playing

  bool running = true;
  SDL_Event event;
  while (running) {
    while (SDL_PollEvent(&event)) {
      switch ( event.type ) {
	case SDL_KEYDOWN:
	  switch (event.key.keysym.sym) {
	  case SDLK_ESCAPE:
	    running = false;
	    break;
	  case SDLK_SPACE:
	    //string->Pluck(0.5);
	    for ( unsigned int k = 0; k < strings.size(); k++ ) {
	      strings[k]->Pluck(0.5);
	    }
	    break;

	  case SDLK_q:
	    strings[0]->Pluck(0.5);
	    break;
	  case SDLK_w:
	    strings[1]->Pluck(0.5);
	    break;
	  case SDLK_e:
	    strings[2]->Pluck(0.5);
	    break;
	  case SDLK_r:
	    strings[3]->Pluck(0.5);
	    break;

	  case SDLK_a:
	    strings[4]->Pluck(0.5);
	    break;
	  case SDLK_s:
	    strings[5]->Pluck(0.5);
	    break;
	  case SDLK_d:
	    strings[6]->Pluck(0.5);
	    break;
	  case SDLK_f:
	    strings[7]->Pluck(0.5);
	    break;

	  case SDLK_z:
	    strings[8]->Pluck(0.5);
	    break;
	  case SDLK_x:
	    strings[9]->Pluck(0.5);
	    break;
	  case SDLK_c:
	    strings[10]->Pluck(0.5);
	    break;
	  case SDLK_v:
	    strings[11]->Pluck(0.5);
	    break;

	  default:
	    break;
	  }
	  break;
      case SDL_QUIT:
	running = false;
	break;
      }
      SDL_Delay(1);
    }
    SDL_Delay(1);
  }
  SDL_Quit();
  return 0;
}
