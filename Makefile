# Top-level build entry point.
#
# Typical invocation:
#   make -j12 CC="ccache gcc" CXX="ccache g++" LGBM=1 LGBM_PATH=/opt/homebrew
#
# User-tunable settings such as ARCH, LGBM, FFTW, ZLIB, PYBIND11_PATH, and
# compiler overrides are defined in or consumed by Makefile.inc.

include Makefile.inc


##############################################################################
# Target Selection
#
# This section controls which binaries/libraries are built by `make all`.
# Safe to change:
# - TARGETS, if you want `all` to build a different default subset.
# Standard/internal:
# - The ARCH / WASM branching reflects platform-specific support choices.
##############################################################################

ifeq ($(ARCH),WINDOWS)
TARGETS = luna destrat behead
else
ifdef WASM
TARGETS = luna libluna destrat
else
TARGETS = luna libluna destrat behead fixrows
endif
endif


##############################################################################
# Source Lists
#
# Standard/internal:
# - These lists define the core source inventory for the main binary and
#   library builds.
# Safe to change:
# - Add or remove source directories here when the codebase layout changes.
##############################################################################

SRCS = globals.cpp eval.cpp param.cpp cmddefs.cpp optdefs.cpp dummy.cpp \
       $(wildcard tests/*.cpp) \
       $(wildcard edf/*.cpp) \
       $(wildcard edfz/*.cpp) \
       $(wildcard defs/*.cpp) \
       $(wildcard tinyxml/*.cpp) \
       $(wildcard helper/*.cpp) \
       $(wildcard timeline/*.cpp) \
       $(wildcard annot/*.cpp) \
       $(wildcard dsp/*.cpp) \
       $(wildcard miscmath/*.cpp) \
       $(wildcard artifacts/*.cpp) \
       $(wildcard spectral/*.cpp) \
       $(wildcard spectral/mtm/*.cpp) \
       $(wildcard spindles/*.cpp) \
       $(wildcard dynamics/*.cpp) \
       $(wildcard intervals/*.cpp) \
       $(wildcard resp/*.cpp) \
       $(wildcard fftw/*.cpp) \
       $(wildcard cwt/*.cpp) \
       $(wildcard stats/*.cpp) \
       $(wildcard stats/catch22/*.cpp) \
       $(wildcard staging/*.cpp) \
       $(wildcard suds/*.cpp) \
       $(wildcard db/*.cpp) \
       $(wildcard ica/*.cpp) \
       $(wildcard clocs/*.cpp) \
       $(wildcard pdc/*.cpp) \
       $(wildcard sstore/*.cpp) \
       $(wildcard dsp/libsamplerate/*.cpp) \
       $(wildcard pops/*.cpp) \
       $(wildcard assoc/*.cpp) \
       $(wildcard lgbm/*.cpp) \
       $(wildcard web/*.cpp) \
       $(wildcard models/*.cpp) \
       $(wildcard lunapi/*.cpp)

CSRCS = $(wildcard db/*.c) \
        $(wildcard stats/*.c) \
        $(wildcard stats/catch22/*.c) \
        $(wildcard dsp/*.c) \
        $(wildcard zlib-1.3.1/*.c) \
        $(wildcard dsp/libsamplerate/*.c)

OBJS = $(SRCS:.cpp=.o) $(CSRCS:.c=.o)
DEPS := $(OBJS:.o=.d)


##############################################################################
# Aggregate Targets
#
# Standard/internal:
# - `all` is the default entry point and expands through TARGETS above.
##############################################################################

all: $(TARGETS)

libluna: libluna.a $(SHARED_LIB)


##############################################################################
# Main Program And Libraries
#
# Standard/internal:
# - These rules describe how the primary executable and library artifacts are
#   linked from the shared object list.
# Safe to change:
# - Link lines only when dependency/link-order requirements actually change.
##############################################################################

luna: main.o $(OBJS)
	$(CXX) -o $@ $^ $(DEP_LIB) $(LDFLAGS)

# Auto-generated header dependency files produced by `-MMD -MP`.
-include $(DEPS)

ifeq ($(ARCH),MAC)
$(SHARED_LIB): $(OBJS)
	$(LD) -dynamiclib $(DEP_LIB) $(LDFLAGS) -o $(SHARED_LIB) $(OBJS)
else ifeq ($(ARCH),LINUX)
$(SHARED_LIB): $(OBJS)
	$(LD) -shared $(DEP_LIB) $(LDFLAGS) -o $(SHARED_LIB) $(OBJS)
endif

libluna.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@


##############################################################################
# Optional / Utility Targets
#
# Safe to change:
# - Add or remove helper utilities as the project evolves.
# Standard/internal:
# - These are intentionally not all part of TARGETS on every platform.
##############################################################################

static: main.o $(OBJS) $(FFTW)/lib/libfftw3.a
	$(CXX) -static -static-libgcc -static-libstdc++ -o luna-static $^

destrat: utils/reader.o libluna.a
	$(CXX) -o $@ $^ -L. $(LDFLAGS) $(DEP_LIB)

regional: utils/region-annotate.o
	$(CXX) -o $@ $^

tocol: utils/tocol.o
	$(CXX) -o $@ $^ $(LDFLAGS)

fixrows: utils/fixrows.o
	$(CXX) -o $@ $^ $(LDFLAGS)

cgi-mapper: utils/cgi-mapper.o libluna.a
	$(CXX) -o $@ $^ $(LDFLAGS) $(DEP_LIB)

behead: utils/behead.o
	$(CXX) -o $@ $^ $(LDFLAGS)

dmerge: utils/merge.o utils/merge-helpers.o
	$(CXX) -o $@ $^

simassoc: utils/simassoc.o libluna.a
	$(CXX) -o $@ $^ $(LDFLAGS) $(DEP_LIB)


##############################################################################
# Test And Housekeeping Targets
#
# Safe to change:
# - The test command arguments, if the project test harness changes.
# Standard/internal:
# - `clean` should continue removing generated build artifacts only.
##############################################################################

.PHONY: all clean libluna static test test-verbose

test: luna
	./luna __LUNA_TESTS__ all

test-verbose: luna
	./luna __LUNA_TESTS__ all verbose

clean:
	-$(RM) $(TARGETS) regional tocol cgi-mapper dmerge simassoc \
		libluna.dylib libluna.so libluna.a luna-static main.o \
		$(OBJS) $(DEPS) $(addsuffix ~,$(SRCS) $(CSRCS))
	-$(RM) utils/*.o utils/*.d utils/*~
