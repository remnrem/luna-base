
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

DEP_LIB += -lfftw3

ifdef LGBM_PATH
DEP_INCLUDES +=-I$(LGBM_PATH)/include
DEP_LIB +=-L$(LGBM_PATH)/lib -L$(LGBM_PATH)
endif





##
## Compiler flags
##

CXX = g++
CC = gcc

CXXFLAGS=-O2 -std=gnu++17 -MMD -MP -I. -I.. -Wno-psabi -Wno-unused-result
# debug
#CXXFLAGS=-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -MMD -MP -I. -I.. -std=gnu++17

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
## Python (lunapi) flags
##

ifdef PYBIND11_PATH
DEP_INCLUDES +=-I$(PYBIND11_PATH)/include
CXXFLAGS+=-DHAS_LUNAPI
endif

ifdef PYTHON_PATH
DEP_INCLUDES += -I$(PYTHON_PATH)
endif

##
## Linker flags
##

LD = g++

LDFLAGS = -L. -L..
# debug
#LDFLAGS=-fsanitize=address,undefined -L. -L..

# assuming mingw64 for windows build:

ifeq ($(ARCH),WINDOWS)
 LDFLAGS += -L/usr/local/lib -L/mingw64/lib/ -L/mingw64/bin/  -static -static-libgcc -static-libstdc++
else
 LDFLAGS += -L/usr/local/lib 
endif

ifdef LGBM_PATH
ifeq ($(ARCH),WINDOWS)
 DEP_LIB += -L/mingw64/bin/
else
 DEP_LIB += -Wl,-rpath,${LGBM_PATH}/lib/ -Wl,-rpath,${LGBM_PATH}/
endif
endif

ifeq ($(LGBM),1)
ifeq ($(ARCH),WINDOWS)
 DEP_LIB += -l:lib_lightgbm.dll 
else
 DEP_LIB += -l_lightgbm 
endif
endif

##
## lomp (needed for static lgbm linking)
##

ifeq ($(LIBOMP),1)
 DEP_LIB += -lomp
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

ifeq ($(ARCH),LINUX)
 SHARED_LIB=libluna.so
endif

ifeq ($(ARCH),MAC)
 SHARED_LIB=libluna.dylib
endif

##
## Misc
##

export CFLAGS
export CXXFLAGS
export CFLAGS
export LD
export LDFLAGS
export LDFLAGS_UTILS

AR = ar
ARFLAGS = rv
RANLIB= ranlib
RM = rm -f
ECHO = echo
SHELL = /bin/sh



