# Top-level build entry point.
#
# Typical invocation:
#   make -j12 CC="ccache gcc" CXX="ccache g++" LGBM=1 LGBM_PATH=/opt/homebrew
#
# User-tunable settings such as ARCH, LGBM, FFTW, ZLIB, PYBIND11_PATH, and
# compiler overrides are defined in or consumed by Makefile.inc.

.DEFAULT_GOAL := all
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
TARGETS = luna destrat destrat-sqlite behead
else
ifdef WASM
TARGETS = luna libluna destrat destrat-sqlite
else
TARGETS = luna libluna destrat destrat-sqlite behead fixrows
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
       $(wildcard dumper/*.cpp) \
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

DUCKDB ?= 0
# The default reader always supports DuckDB; DUCKDB controls Luna's writer.
DUCKDB_PATH ?= /usr/local
DUCK_INCLUDES = -I"$(DUCKDB_PATH)" -I"$(DUCKDB_PATH)/include"
ifeq ($(ARCH),WINDOWS)
DUCK_LIB_DIRS = $(wildcard $(DUCKDB_PATH) $(DUCKDB_PATH)/bin)
DUCK_LIBS = $(foreach dir,$(DUCK_LIB_DIRS),-L"$(dir)") -l:duckdb.dll
else
DUCK_LIB_DIRS = $(wildcard $(DUCKDB_PATH) $(DUCKDB_PATH)/lib)
DUCK_LIBS = $(foreach dir,$(DUCK_LIB_DIRS),-L"$(dir)" -Wl,-rpath,"$(dir)") -lduckdb
endif
duckdb/%.o: CPPFLAGS += $(DUCK_INCLUDES)
ifeq ($(DUCKDB),1)
DUCK_WRITER = duckdb/writer.o duckdb/catalog.o
main.o: CPPFLAGS += -DHAS_DUCKDB $(DUCK_INCLUDES)
endif

all: $(TARGETS)

# Only main.o changes with the feature flag. Force its rebuild on a mode switch.
DUCK_MODE_STAMP = duckdb/.build-mode-$(DUCKDB)
$(DUCK_MODE_STAMP):
	$(RM) duckdb/.build-mode-*
	touch $@
main.o: $(DUCK_MODE_STAMP)

# Keep the old DuckDB reader name available for existing scripts.
destrat destrat-duckdb: duckdb/reader.o duckdb/reader_source.o duckdb/reader_sqlite.o duckdb/reader_typed.o duckdb/catalog.o libluna.a
	$(CXX) -o $@ $^ $(DEP_LIB) $(LDFLAGS) $(DUCK_LIBS)

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

luna: main.o $(OBJS) $(DUCK_WRITER)
	$(CXX) -o $@ $^ $(DEP_LIB) $(LDFLAGS) $(DUCK_LIBS)

# Auto-generated header dependency files produced by `-MMD -MP`.
-include main.d $(DEPS)
-include utils/reader2.d duckdb/writer.d duckdb/reader.d duckdb/reader_source.d duckdb/reader_sqlite.d duckdb/reader_typed.d duckdb/catalog.d

.PHONY: hdstats-pc-schema

hdstats-pc-schema:
	python3 scripts/generate_hdstats_pc_schema.py docs/hdstats/hdstats-pc.spec.json pops/hdstats-pc-schema.h

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

# The original utils/reader.cpp is retained as reference source only.
destrat-sqlite: utils/reader2.o libluna.a
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

.PHONY: all clean libluna static test test-verbose test-destrat test-destrat2

test: luna
	./luna __LUNA_TESTS__ all

test-verbose: luna
	./luna __LUNA_TESTS__ all verbose

test-destrat: destrat
	python3 tests/test_destrat2.py

# Compatibility alias for the existing test command; builds only destrat.
test-destrat2: test-destrat

clean:
	-$(RM) $(TARGETS) destrat2 regional tocol cgi-mapper dmerge simassoc \
		libluna.dylib libluna.so libluna.a luna-static main.o main.d \
		$(OBJS) $(DEPS) $(addsuffix ~,$(SRCS) $(CSRCS))
	-$(RM) utils/*.o utils/*.d utils/*~ duckdb/*.o duckdb/*.d duckdb/.build-mode-* destrat-duckdb
