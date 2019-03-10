
##
## ARCH is LINUX, MAC or WINDOWS
##

ARCH=LINUX

##
## If dependencies (FFTW, libz if needed, and optionally, libhpdf and libpng) are installed locally, 
## the DEP variable will be set on the command line
##

ifdef LUNA_DEP
DEP_INCLUDES = -I$(LUNA_DEP)/include
DEP_LIB = -L$(LUNA_DEP)/lib
endif


##
## Basic C/C++ compilation flags
## 

CXX = g++

CXXFLAGS = -O2 -I. -I.. -I/usr/local/include $(DEP_INCLUDES) -Wno-logical-op-parentheses -Wno-deprecated-register

CFLAGS = -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION -DT_LINUX

##
## By default, do not try to add pdf_t support (i.e. omit -lhpdf and -lpng)
##

CXXFLAGS += -DNO_HPDFLIB

##
## Static binary: 'make static STATIC=1', otherwise set PIC code
##

ifndef STATIC
CXXFLAGS += -fPIC
CFLAGS += -fPIC
endif


##
## Linker flags
##

LD = g++ 
LDFLAGS = $(DEP_LIB) -L/usr/local/lib


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
