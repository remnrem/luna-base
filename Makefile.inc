
# ARCH is LINUX or MAC 
ARCH=LINUX

# if dependencies are installed locally, the DEP variable will be set on the command line
ifdef LUNA_DEP
DEP_INCLUDES = -I$(LUNA_DEP)/include
DEP_LIB = -L$(LUNA_DEP)/lib
endif

CXX = g++
CXXFLAGS = -O2 -I. -I.. -I/usr/local/include $(DEP_INCLUDES) -std=c++0x -Wno-logical-op-parentheses -Wno-deprecated-register -fPIC
CFLAGS = -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -DT_LINUX -fPIC

LD = g++ 
LDFLAGS = $(DEP_LIB) -L/usr/local/lib 

export CFLAGS
export CXXFLAGS
export CFLAGS
export LD
export LDFLAGS

AR = ar
ARFLAGS =
RANLIB= ranlib
RM = rm -f
ECHO = echo
SHELL = /bin/sh

.SILENT :
