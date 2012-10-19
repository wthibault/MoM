HOST_PLATFORM := $(shell $(CPP) -dumpmachine)
$(info $(HOST_PLATFORM))

ifeq   "$(HOST_PLATFORM)" "i686-apple-darwin11"
RTAUDIO_DIR = /Users/tebo/Software/rtaudio-4.0.10
CCFLAGS = -fopenmp -O3
endif

ifeq   "$(HOST_PLATFORM)" "x86_64-linux-gnu"
#RTAUDIO_DIR = /home/tebo/Software/rtaudio-4.0.10
CCFLAGS = -fopenmp -O3
endif

all : rtstring viz string string-omp

rtstring : rtstring.cpp
	g++ $(CCFLAGS) rtstring.cpp -o rtstring \
		-I$(RTAUDIO_DIR) -L$(RTAUDIO_DIR) -lrtaudio \
		-lpthread

viz : viz.cpp util.o
	g++ $(CCFLAGS) viz.cpp -o viz \
		util.o \
		-I$(RTAUDIO_DIR) -L$(RTAUDIO_DIR) -lrtaudio \
		-lglut -lGL -lGLU \
		-lpthread

util.o : util.cpp

string : string.cpp
	g++ $(CCFLAGS) string.cpp -o string -L$(RTAUDIO_DIR) -lrtaudio

string-omp: string-omp.cpp
	g++ $(CCFLAGS) string-omp.cpp -o string-omp -L$(RTAUDIO_DIR) -lrtaudio

clean :
	rm -f a.out
	rm -f *.o
	rm -f *~
	rm -f string
	rm -f viz
	rm -f rtstring
	rm -f string-omp
	rm -f *.raw
