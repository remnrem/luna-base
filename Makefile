include Makefile.inc
ifeq ($(ARCH),WINDOWS)
  TARGETS = luna destrat behead
else
  ifdef WASM
    TARGETS=luna libluna destrat
  else
    TARGETS = luna libluna destrat behead fixrows 
  endif
endif

SRCS = globals.cpp eval.cpp param.cpp cmddefs.cpp optdefs.cpp dummy.cpp \
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
        $(wildcard dsp/*.c) \
        $(wildcard zlib-1.3/*.c) \
        $(wildcard dsp/libsamplerate/*.c)

OBJS = $(SRCS:.cpp=.o) $(CSRCS:.c=.o) 

DEPS := $(OBJS:.o=.d)

#
# targets
#

all: $(TARGETS)

luna: main.o $(OBJS)
	$(CXX) -o $@ $^ $(DEP_LIB) $(LDFLAGS)

libluna: libluna.a $(SHARED_LIB)

# header dependencies

-include $(DEPS)

# shared library (libluna)

ifeq ($(ARCH),MAC)
$(SHARED_LIB) : $(OBJS)
	$(LD) -dynamiclib $(DEP_LIB) $(LDFLAGS) -o $(SHARED_LIB) $(OBJS)
else ifeq ($(ARCH),LINUX)
$(SHARED_LIB) : $(OBJS)
	$(LD) -shared $(DEP_LIB) $(LDFLAGS) -o $(SHARED_LIB) $(OBJS)
endif

# objects
libluna.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

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
	$(CXX)  -o $@ $^ $(LDFLAGS) $(DEP_LIB)

behead: utils/behead.o
	$(CXX) -o $@ $^  $(LDFLAGS)

dmerge: utils/merge.o utils/merge-helpers.o
	$(CXX) -o $@ $^ 

simassoc: utils/simassoc.o libluna.a
	$(CXX) -o $@ $^  $(LDFLAGS) $(DEP_LIB)

.PHONY: clean

clean:
	-$(RM) $(TARGETS) libluna.dylib libluna.so libluna.a main.o $(OBJS) $(DEPS) $(addsuffix ~,$(SRCS) $(CSRCS))
	-$(RM) utils/*.o utils/*.d utils/*~
