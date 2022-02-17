
##
## ARCH is LINUX, MAC or WINDOWS
##

ifndef ARCH
ARCH=LINUX
endif

##
## If dependencies (FFTW, LGBM) are installed locally, 
## the DEP variable will be set on the command line
##

DEP_INCLUDES=
DEP_LIB=

ifdef ZLIB
DEP_INCLUDES += -I$(ZLIB)/include
DEP_LIB += -L$(ZLIB)/lib
endif

ifdef FFTW
DEP_INCLUDES += -I$(FFTW)/include
DEP_LIB = -L$(FFTW)/lib
ifeq ($(ARCH),LINUX)
DEP_LIB	+= -Wl,-rpath,$(FFTW)/lib
endif
endif

ifdef LGBM_PATH
DEP_INCLUDES +=-I$(LGBM_PATH)/include
DEP_LIB =-L$(LGBM_PATH)
endif


##
## Compiler flags
##

CXX = g++
CC = gcc

CXXFLAGS=-O2 -std=gnu++11 -MMD -MP -I. -I..
#CXXFLAGS=-O0 -std=gnu++11 -MMD -MP -I. -I.. -pg

CXXFLAGS +=$(DEP_INCLUDES)

ifeq ($(ARCH),MAC)
CXXFLAGS+=-I/usr/local/include
endif

ifeq ($(ARCH),WINDOWS)
CXXFLAGS+=-DWINDOWS
endif

ifeq ($(LGBM),1)
CXXFLAGS+=-DHAS_LGBM
endif

CFLAGS=-O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION


##
## Linker flags
##

LD = g++ 

LDFLAGS = $(DEP_LIB) -L. -L..

# assuming mingw64 for windows build:

ifeq ($(ARCH),WINDOWS)
LDFLAGS += -L/usr/local/lib -L/mingw64/lib/ -L/mingw64/bin/ -lfftw3 -lz -static -static-libgcc -static-libstdc++
else
LDFLAGS += -L/usr/local/lib -lfftw3 -lz
endif

ifeq ($(LGBM),1)
ifeq ($(ARCH),LINUX)
LDFLAGS += -Wl,-rpath=${LGBM_PATH}
endif
LDFLAGS += -l_lightgbm
endif



##
## Static binary: 'make static STATIC=1', otherwise set PIC code
##

ifdef STATIC
LDFLAGS += -static
ifdef LINUX
LDFLAGS += -static-libgcc -static-libstdc++
endif
endif

##
## Build with dynamic linking, and libluna.so (default)
##

ifndef STATIC
CXXFLAGS+=-fPIC
CFLAGS+=-fPIC
LDFLAGS+=-fPIC
endif

##
## Shared library name
##

SHARED_LIB=libluna.dylib

##
## Misc
##

export CFLAGS
export CXXFLAGS
export CFLAGS
export LD
export LDFLAGS

AR = ar
ARFLAGS = rv
RANLIB= ranlib
RM = rm -f
ECHO = echo
SHELL = /bin/sh
