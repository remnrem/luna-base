##############################################################################
# Platform Selection
#
# Safe to change:
# - ARCH on the command line or in the environment, e.g. `make ARCH=MAC`.
# Standard/internal:
# - Supported values are LINUX, MAC, and WINDOWS.
##############################################################################

ifndef ARCH
ARCH = LINUX
endif


##############################################################################
# User-Configurable Build Options
#
# Safe to change:
# - These are the main knobs expected to be overridden by users or CI:
#   `make LGBM=0`, `make LGBM_PATH=/opt/homebrew`, `make STATIC=1`, etc.
# Standard/internal:
# - Defaults here should stay conservative and portable.
##############################################################################

LGBM ?= 1
LGBM_PATH ?= /opt/homebrew


##############################################################################
# External Dependency Discovery
#
# Safe to change:
# - FFTW, ZLIB, PYBIND11_PATH, and PYTHON_PATH when dependencies live outside
#   standard system locations.
# Standard/internal:
# - DEP_INCLUDES / DEP_LIB are internal accumulators used by later flag blocks.
##############################################################################

DEP_INCLUDES =
DEP_LIB =

ifdef ZLIB
DEP_INCLUDES += -I$(ZLIB)/include
DEP_LIB += -L$(ZLIB)/lib
endif

ifdef FFTW
DEP_INCLUDES += -I$(FFTW)/include
DEP_LIB += -L$(FFTW)/lib
ifeq ($(ARCH),LINUX)
DEP_LIB += -Wl,-rpath,$(FFTW)/lib
endif
endif

# FFTW is treated as a standard dependency for normal builds.
DEP_LIB += -lfftw3

ifeq ($(LGBM),1)
ifdef LGBM_PATH
DEP_INCLUDES += -I$(LGBM_PATH)/include
DEP_LIB += -L$(LGBM_PATH)/lib -L$(LGBM_PATH)
endif
endif


##############################################################################
# Compiler Toolchain And Compile Flags
#
# Safe to change:
# - CC, CXX, optimization/debug flags, warning flags, and feature defines.
# Standard/internal:
# - `-MMD -MP` is required for the auto-generated dependency files included by
#   the top-level Makefile.
##############################################################################

CXX = g++
CC = gcc

CXXFLAGS = -O2 -std=gnu++17 -MMD -MP -I. -I.. -Wno-psabi -Wno-unused-result
# Debug alternative:
# CXXFLAGS = -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
#            -MMD -MP -I. -I.. -std=gnu++17

ifeq ($(ARCH),MAC)
CXXFLAGS += -I/usr/local/include
endif

ifeq ($(ARCH),WINDOWS)
CXXFLAGS += -DWINDOWS
endif

ifeq ($(LGBM),1)
CXXFLAGS += -DHAS_LGBM
endif

CFLAGS = -O2 -DSQLITE_THREADSAFE=0 -DSQLITE_OMIT_LOAD_EXTENSION


##############################################################################
# Python / lunapi Support
#
# Safe to change:
# - PYBIND11_PATH and PYTHON_PATH to point at the desired Python headers.
# Standard/internal:
# - HAS_LUNAPI is enabled only when pybind11 headers are explicitly supplied.
##############################################################################

ifdef PYBIND11_PATH
DEP_INCLUDES += -I$(PYBIND11_PATH)/include
CXXFLAGS += -DHAS_LUNAPI
endif

ifdef PYTHON_PATH
DEP_INCLUDES += -I$(PYTHON_PATH)
endif

CXXFLAGS += $(DEP_INCLUDES)


##############################################################################
# Linker Selection And Link Flags
#
# Safe to change:
# - LDFLAGS, LD, and library search paths when toolchain layout differs.
# Standard/internal:
# - DEP_LIB is extended here with runtime search paths and optional libraries.
##############################################################################

LD = g++

LDFLAGS = -L. -L..
# Debug alternative:
# LDFLAGS = -fsanitize=address,undefined -L. -L..

# Windows assumptions currently target a mingw64-style environment.
ifeq ($(ARCH),WINDOWS)
LDFLAGS += -L/usr/local/lib -L/mingw64/lib/ -L/mingw64/bin/ -static \
           -static-libgcc -static-libstdc++
else
LDFLAGS += -L/usr/local/lib
endif

ifeq ($(LGBM),1)
ifdef LGBM_PATH
ifeq ($(ARCH),WINDOWS)
DEP_LIB += -L/mingw64/bin/
else
DEP_LIB += -Wl,-rpath,$(LGBM_PATH)/lib/ -Wl,-rpath,$(LGBM_PATH)/
endif
endif
endif

ifeq ($(LGBM),1)
ifeq ($(ARCH),WINDOWS)
DEP_LIB += -l:lib_lightgbm.dll
else
DEP_LIB += -l_lightgbm
endif
endif


##############################################################################
# Optional OpenMP Support
#
# Safe to change:
# - LIBOMP=1 when static LightGBM builds require an explicit OpenMP runtime.
##############################################################################

ifeq ($(LIBOMP),1)
DEP_LIB += -lomp
endif


##############################################################################
# Static vs Dynamic Build Mode
#
# Safe to change:
# - STATIC=1 to request a static binary-oriented build.
# Standard/internal:
# - Default builds are PIC-enabled so `libluna.so` / `libluna.dylib` can be
#   produced from the shared object list.
##############################################################################

ifdef STATIC
LDFLAGS += -static
ifeq ($(ARCH),LINUX)
LDFLAGS += -static-libgcc -static-libstdc++
endif
endif

ifndef STATIC
CXXFLAGS += -fPIC
CFLAGS += -fPIC
LDFLAGS += -fPIC
endif


##############################################################################
# Shared Library Naming
#
# Standard/internal:
# - The top-level Makefile uses SHARED_LIB when building `libluna`.
##############################################################################

ifeq ($(ARCH),LINUX)
SHARED_LIB = libluna.so
endif

ifeq ($(ARCH),MAC)
SHARED_LIB = libluna.dylib
endif


##############################################################################
# Exported Tools And Helpers
#
# Standard/internal:
# - These exports/variables are consumed by sub-makes and the top-level file.
# Safe to change:
# - Tool locations such as AR, RANLIB, RM, or SHELL when needed.
##############################################################################

export CFLAGS
export CXXFLAGS
export LD
export LDFLAGS

AR = ar
ARFLAGS = rv
RANLIB = ranlib
RM = rm -f
ECHO = echo
SHELL = /bin/sh
