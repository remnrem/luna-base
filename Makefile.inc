
##
## ARCH is LINUX, MAC or WINDOWS
##

ARCH=LINUX


##
## Basic C/C++ compilation flags
## 

CXX = g++

CXXFLAGS = -O2 -I. -I.. -Wno-logical-op-parentheses -Wno-deprecated-register

CXXFLAGS += $(DEP_INCLUDES) 

ifdef MAC
CXXFLAGS += -I/usr/local/include
endif

CFLAGS = -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION


##
## If dependencies (FFTW, if needed, and optionally, libhpdf and libpng) are installed locally, 
## the DEP variable will be set on the command line
##

DEP_INCLUDES = 
DEP_LIB = 


ifdef FFTW
DEP_INCLUDES += -I$(FFTW)/include
DEP_LIB = -L$(FFTW)/lib
endif

ifdef HPDF
DEP_INCLUDES += -I$(HPDF)/include
DEP_LIB = -L$(HPDF)/lib
else
CXXFLAGS += -DNO_HPDFLIB 
endif

ifdef LIBPNG
DEP_INCLUDES += -I$(LIBPNG)/include
DEP_LIB = -L$(LIBPNG)/lib
else
CXXFLAGS += -DNO_HPDFLIB 
endif



##
## Linker flags
##

LD = g++ 
LDFLAGS = $(DEP_LIB)

ifndef WINDOWS
LDFLAGS += -L/usr/local/lib
endif



##
## Static binary: 'make static STATIC=1', otherwise set PIC code
##

ifndef LUNA_STATIC
CXXFLAGS += -fPIC
CFLAGS += -fPIC
else
LDFLAGS += -static -static-libgcc -static-libstdc++
endif


##
## Misc
##

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
