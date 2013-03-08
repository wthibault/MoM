#OPTS = -g
OPTS = -O3 -DUSE_HISTOGRAMS -DREFRESH_PERIOD=20

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
OPTS +=  -DNUM_MASSES=500
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


all : 	string-apparatus

string-apparatus : string-apparatus.cpp \
		StringModel.o \
		ui.o \
		$(SSG_INCS) \
		$(LIB_DIR)/ssg.a \
		$(SSG_DATA)/objfiles/string-apparatus.obj 
	g++ $(OPTS) -o string-apparatus $(INCLUDES) string-apparatus.cpp \
		StringModel.o \
		ui.o \
		$(LIB_DIR)/ssg.a \
		-I$(RTAUDIO_DIR)  \
		$(LIBS)


StringModel.o : StringModel.h StringModel.cpp
	g++ $(OPTS) -c StringModel.cpp $(INCLUDES)

ui.o : ui.h ui.cpp 
	g++ $(OPTS) -c ui.cpp $(INCLUDES)

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
