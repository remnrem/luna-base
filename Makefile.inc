
CXX = g++
#CXX = clang++
CXXFLAGS = -O2 -I.. -I/usr/local/include -std=c++0x -Wno-logical-op-parentheses -Wno-deprecated-register 
CFLAGS = -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -DT_LINUX

LD = g++ 
LDLAGS = -L/usr/local/lib

# if dependencies are installed locally, denote here
CXXFLAGS += -I/Users/shaun/src/luna-depends/depends/include
LDFLAGS += -L/Users/shaun/src/luna-depends/depends/lib

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
