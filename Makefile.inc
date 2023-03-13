
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
DEP_LIB +=-L$(LGBM_PATH)
endif


##
## Compiler flags
##

ifeq ($(WASM),1)
 CXX = em++
 CC = emcc
else
 CXX = g++
 CC = gcc
endif


CXXFLAGS=-O2 -std=gnu++11 -MMD -MP -I. -I.. 
#CXXFLAGS=-O0 -std=gnu++11 -MMD -MP -I. -I.. -pg

CXXFLAGS +=$(DEP_INCLUDES)

ifeq ($(WASM),1)
 CXXFLAGS+=-DWASM -s USE_ZLIB=1
endif


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

ifeq ($(WASM),1)
 LD = em++
else
 LD = g++
endif


LDFLAGS = $(DEP_LIB) -L. -L..

# assuming mingw64 for windows build:

ifeq ($(ARCH),WINDOWS)
 LDFLAGS += -L/usr/local/lib -L/mingw64/lib/ -L/mingw64/bin/ -lfftw3 -lz -static -static-libgcc -static-libstdc++
else
 ifeq ($(WASM),1)
  LDFLAGS += -L${WASM_DEPS} -lfftw3 -s USE_ZLIB=1
 else
  LDFLAGS += -L/usr/local/lib -lfftw3 -lz
 endif
endif

ifdef LGBM_PATH
ifeq ($(ARCH),WINDOWS)
LDFLAGS += -L/mingw64/bin/
else
LDFLAGS += -Wl,-rpath,${LGBM_PATH}
endif
endif

ifeq ($(LGBM),1)
ifeq ($(ARCH),WINDOWS)
LDFLAGS += -l:lib_lightgbm.dll
else
LDFLAGS += -l_lightgbm
endif
endif


##
## WASM flags
##

ifdef WASM
#LDFLAGS += -s LINKABLE=1 -s EXPORT_ALL=1 -s ERROR_ON_UNDEFINED_SYMBOLS=0 -sALLOW_MEMORY_GROWTH -sFORCE_FILESYSTEM  --preload-file tutorial/
 LDFLAGS += -s LINKABLE=1 -s EXPORT_ALL=1 -s ERROR_ON_UNDEFINED_SYMBOLS=0 -sALLOW_MEMORY_GROWTH -sFORCE_FILESYSTEM
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


ifdef WASM
 AR = emar
 ARFLAGS = rv
 RANLIB= emranlib
else
 AR = ar
 ARFLAGS = rv
 RANLIB= ranlib
endif

RM = rm -f
ECHO = echo
SHELL = /bin/sh
