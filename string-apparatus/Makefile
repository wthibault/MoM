#OPTS = -g
OPTS = -O3 -DLOCK_STRING
#OPTS = -O3

# change these for your system
INC_DIR=../ssg/include
LIB_DIR=../ssg/src
SSG_DATA=../ssg/data

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
	-L/usr/X11/lib \
	/usr/local/lib/libftgl.a \
	-L/usr/local/lib -lfftw3 \
	-headerpad_max_install_names
OPTS +=  -DNUM_MASSES=256
POST_BUILD_STEP = cp string-apparatus macosx/StringApparatus.app/Contents/MacOS
endif


# for MINGW on windows
ifeq   "$(HOST_PLATFORM)" "mingw32"
RTAUDIO_DIR = rtaudio-4.0.10
INCLUDES += -I$(INC_DIR) \
	-I$(RTAUDIO_DIR) \
	-I../include \
	-I/usr/local/include \
	-I../ssg/win32/ssg-lib/glew-1.9.0/include  -DGLEW_STATIC  \
	-I../ssg/win32/ssg-lib/glut-3.7.6-bin \
	-D__WIN32__


LIBS = $(RTAUDIO_DIR)/librtaudio.a \
	/usr/local/lib/libfftw3.a \
	/usr/local/lib/libfltk.a \
	-ldsound -lwinmm\
	-L../ssg/win32/ssg-lib/glut-3.7.6-bin -lglut32

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
all : 	string-apparatus 

string-apparatus : main.cpp \
		StringModelPrimitive.h \
		FFTPrimitive.h \
		StringModel.o \
		gui.o \
		$(SSG_INCS) \
		$(LIB_DIR)/ssg.a 
#		$(SSG_DATA)/objfiles/string-apparatus.obj 
	g++ $(OPTS) -o string-apparatus $(INCLUDES) \
		main.cpp \
		StringModel.o \
		gui.o \
		$(LIB_DIR)/ssg.a \
		-I$(RTAUDIO_DIR)  \
		$(LIBS) \
		$(LINKFLTK_ALL)
	$(POST_BUILD_STEP)


StringModel.o : StringModel.h StringModel.cpp
	g++ $(OPTS) -c StringModel.cpp $(INCLUDES)

gui.o : gui.h gui.cpp 
	g++ $(OPTS) -c gui.cpp $(INCLUDES)

objfiles/string-apparatus.obj : blender/string-apparatus.obj
	cp blender/string-apparatus.obj objfiles
	cp blender/string-apparatus.mtl objfiles

objfiles/string-scene.obj : blender/string-scene.obj
	cp blender/string-scene.obj objfiles
	cp blender/string-scene.mtl objfiles

ifeq   "$(HOST_PLATFORM)" "mingw32"
package :
	rm -rf win32/StringApparatus
	mkdir win32/StringApparatus
	cp string-apparatus.exe win32/StringApparatus
	cp -r objfiles/ win32/StringApparatus/objfiles
	cp -r shaders120/ win32/StringApparatus/shaders120
endif

clean :
	rm -f *.o
	rm -f shaders/*~
	rm -f *~
	rm -f string-apparatus
	rm -rf *.dSYM
	rm -rf win32/StringApparatus
