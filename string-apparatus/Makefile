#OPTS = -g
OPTS = -O3

# change these for your system
INC_DIR=../../4848/SSG/include
LIB_DIR=../../4848/SSG/src
SSG_DATA=../../4848/SSG/data

HOST_PLATFORM := $(shell $(CPP) -dumpmachine)
$(info $(HOST_PLATFORM))

# settings for turing
ifeq   "$(HOST_PLATFORM)" "x86_64-redhat-linux"
INCLUDES = -I/opt/viz/include
LIBS = -lGL -lGLU -lglut -L/opt/viz/lib -lGLEW
endif

ifeq   "$(HOST_PLATFORM)" "x86_64-linux-gnu"
#for sn336,sn148, ubuntu 
RTAUDIO_DIR = /Users/tebo/Software/rtaudio-4.0.10
INCLUDES = -I$(INC_DIR)
LIBS = -lGL -lGLU -lglut -lGLEW -lpthread
endif

# for macosx
# NOTE for mac you have to use the GLSL 1.2 shaders....
ifeq   "$(HOST_PLATFORM)" "i686-apple-darwin11"
RTAUDIO_DIR = /Users/tebo/Software/rtaudio-4.0.10
INCLUDES = -I$(INC_DIR) -framework OpenGL -framework GLUT -DUSE_APPLE_VERTEXARRAYS \
	-I$(RTAUDIO_DIR)
LIBS = -L$(RTAUDIO_DIR)
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
	$(SSG_INCS) \
	$(LIB_DIR)/ssg.a \
	$(SSG_DATA)/objfiles/string-apparatus.obj \
	StringModel.o
	g++ $(OPTS) -o string-apparatus $(INCLUDES) string-apparatus.cpp $(LIB_DIR)/ssg.a $(LIBS) StringModel.o \
		-I$(RTAUDIO_DIR) -L$(RTAUDIO_DIR) -lrtaudio \


StringModel.o : StringModel.h StringModel.cpp
	g++ $(OPTS) -c StringModel.cpp $(INCLUDES)

$(SSG_DATA)/objfiles/string-apparatus.obj : blender/string-apparatus.obj
	cp blender/string-apparatus.obj $(SSG_DATA)/objfiles
	cp blender/string-apparatus.mtl $(SSG_DATA)/objfiles

clean :
	rm -f *.o
	rm -f *~
	rm -f string-apparatus
	rm -rf *.dSYM
