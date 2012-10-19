# for macosx
#FLAGS = -lm -framework OpenGL -framework AGL  -framework IOKit -framework Cocoa -framework QuickTime -framework ApplicationServices -framework Carbon -framework AudioToolbox -framework AudioUnit

LIBS = -L/usr/local/lib -lSDL -lSDLmain

all: pluckfile sdlpluck

sdlpluck : SDLpluck.cpp
	g++ $(FLAGS) -o sdlpluck SDLpluck.cpp $(LIBS)

pluckfile : pluckfile.cpp
	g++ $(FLAGS) -o pluckfile pluckfile.cpp $(LIBS)

clean:
	rm -f sdlpluck
	rm -f pluckfile
	rm -f *~
