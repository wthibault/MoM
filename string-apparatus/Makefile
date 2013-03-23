#OPTS = -g
OPTS = -O3 -DLOCK_STRING
#OPTS = -O3

# change these for your system
INC_DIR=../../4848/ssg/include
LIB_DIR=../../4848/ssg/src
SSG_DATA=../../4848/ssg/data

HOST_PLATFORM := $(shell $(CPP) -dumpmachine)
$(info $(HOST_PLATFORM))

# settings for turing
ifeq   "$(HOST_PLATFORM)" "x86_64-redhat-linux"
INCLUDES = -I/opt/viz/include
LIBS = -lGL -lGLU -lglut -L/opt/viz/lib -lGLEW
endif

ifeq   "$(HOST_PLATFORM)" "x86_64-linux-gnu"
#for sn336,sn148, ubuntu 
# for package libfft3w-dev
#FFTW=/usr/lib/libfftw3.a
# for build from source, with default prefix (/usr/local)
#FFTW=/usr/local/lib/libfftw3.a
# for both , only not static linked...
FFTW = $(shell pkg-config fftw3 --libs)
RTAUDIO_DIR = ./rtaudio-4.0.10
INCLUDES = -I$(INC_DIR) -I$(RTAUDIO_DIR) $(shell pkg-config ftgl --cflags) $(shell pkg-config fftw3 --cflags)
LIBS = -lGL -lGLU -lglut -lGLEW \
	$(shell pkg-config ftgl --libs) \
	-lpthread $(FFTW) \
	$(RTAUDIO_DIR)/librtaudio.a -lasound -lpthread
OPTS +=  -DNUM_MASSES=1000
endif

# for macosx
# NOTE for mac you have to use the GLSL 1.2 shaders....
ifeq   "$(HOST_PLATFORM)" "i686-apple-darwin11"
RTAUDIO_DIR = rtaudio-4.0.10
INCLUDES = -I$(INC_DIR) \
	-framework OpenGL \
	-framework GLUT -DUSE_APPLE_VERTEXARRAYS \
	-I$(RTAUDIO_DIR) \
	-I/usr/X11/include -I/usr/X11/include/freetype2

LIBS = $(RTAUDIO_DIR)/librtaudio.a \
	-lpthread -framework CoreAudio -framework CoreFoundation \
	-L/usr/X11/lib -lftgl \
	-L/usr/local/lib -lfftw3
OPTS +=  -DNUM_MASSES=256
endif

SSG_INCS = $(INC_DIR)/ssg.h \
	$(INC_DIR)/CheckError.h \
	$(INC_DIR)/InitShader.h \
	$(INC_DIR)/Texture.h \
	$(INC_DIR)/TextureLoader.h \
	$(INC_DIR)/ObjLoader.h \
	$(INC_DIR)/FrameBufferObject.h \
	$(INC_DIR)/ParticleSystem.h \
	$(INC_DIR)/Camera.h \
	$(INC_DIR)/Trackball.h 


###### 
FLTKCONFIG   = /usr/local/bin/fltk-config
CXX          = $(shell $(FLTKCONFIG) --cxx)
#CXXFLAGS     = $(shell $(FLTKCONFIG) --cxxflags) -g -Wall -I. -I$(SSG_DIR)/include
INCLUDES    += $(shell $(FLTKCONFIG) --cxxflags)
LINKFLTK     = $(shell $(FLTKCONFIG) --ldstaticflags)
LINKFLTK_GL  = $(shell $(FLTKCONFIG) --use-gl --ldstaticflags)
LINKFLTK_IMG = $(shell $(FLTKCONFIG) --use-images --ldstaticflags)
LINKFLTK_ALL = $(shell $(FLTKCONFIG) --use-images --use-gl --ldstaticflags) 
ifeq "$(MACHINE)" "i686-apple-darwin11"
LINKFLTK_ALL = $(shell $(FLTKCONFIG) --use-images --use-gl --ldstaticflags) \
	       $(SSG_DIR)/src/ssg.a -framework OpenGL
endif

#.SUFFIXES: .cxx .h .fl .o $(EXEEXT)

# HOW TO COMPILE
#.cxx.o:
#	@echo "*** Compile $<..."
#	$(CXX) -I.. $(CXXFLAGS) -c $< -o $@

## HOW TO LINK
#.o$(EXEEXT):
#	@echo "*** Linking $@..."
#	$(CXX) $< $(LINKFLTK_ALL) -o $@

##############
all : 	string-apparatus $(SSG_DATA)/objfiles/string-apparatus.obj $(SSG_DATA)/objfiles/string-scene.obj

string-apparatus : main.cpp \
		StringModelPrimitive.h \
		FFTPrimitive.h \
		StringModel.o \
		gui.o \
		$(SSG_INCS) \
		$(LIB_DIR)/ssg.a \
		$(SSG_DATA)/objfiles/string-apparatus.obj 
	g++ $(OPTS) -o string-apparatus $(INCLUDES) \
		main.cpp \
		StringModel.o \
		gui.o \
		$(LIB_DIR)/ssg.a \
		-I$(RTAUDIO_DIR)  \
		$(LIBS) \
		$(LINKFLTK_ALL)


StringModel.o : StringModel.h StringModel.cpp
	g++ $(OPTS) -c StringModel.cpp $(INCLUDES)

gui.o : gui.h gui.cpp 
	g++ $(OPTS) -c gui.cpp $(INCLUDES)

$(SSG_DATA)/objfiles/string-apparatus.obj : blender/string-apparatus.obj
	cp blender/string-apparatus.obj $(SSG_DATA)/objfiles
	cp blender/string-apparatus.mtl $(SSG_DATA)/objfiles

$(SSG_DATA)/objfiles/string-scene.obj : blender/string-scene.obj
	cp blender/string-scene.obj $(SSG_DATA)/objfiles
	cp blender/string-scene.mtl $(SSG_DATA)/objfiles

clean :
	rm -f *.o
	rm -f shaders/*~
	rm -f *~
	rm -f string-apparatus
	rm -rf *.dSYM
