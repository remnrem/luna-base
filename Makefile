include Makefile.inc

DIRS = edf tinyxml helper timeline annot dsp miscmath spindles	\
artifacts intervals fftw cwt defs zfile stats graphics staging 	\
db ica clocs pdc sstore

EXE	= luna
LUNALIB = libluna.so
OBJS	= main.o globals.o eval.o

OBJLIBS = libdefs.a libedf.a libtinyxml.a libhelper.a libtimeline.a	\
libannot.a libdsp.a libmiscmath.a libspindles.a libartifacts.a		\
libintervals.a libfftwrap.a libcwt.a libzfile.a libstats.a		\
libgraphics.a libstaging.a libdb.a libica.a libclocs.a libpdc.a		\
libsstore.a

LIBS = -L. -lspindles -lica -lannot -ldefs -lartifacts -ledf -lhelper \
 -ltimeline -lstaging -lfftwrap -ldsp -lmiscmath -lintervals \
 -ltinyxml -lcwt -lclocs -lpdc -lzfile -lstats -lgraphics -ldb -lsstore \
 -lfftw3 -lhpdf -lpng -lsamplerate -lz

all : $(EXE) $(LUNALIB) utils

$(EXE) : main.o globals.o eval.o $(OBJLIBS)
	$(ECHO) $(LD) $(LDFLAGS) -o $(EXE) $(OBJS) $(LIBS)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(LUNALIB) : globals.o eval.o $(OBJLIBS)

ifeq ($(ARCH),MAC)
	$(ECHO) "building libluna.dylib..."
	$(LD) -dynamiclib -fPIC $(LDFLAGS) -o libluna.dylib eval.o globals.o  *.a -lz -lfftw3 -lhpdf -lsamplerate
else
	$(ECHO) "building libluna.so..."
	$(LD) -shared     -fPIC $(LDFLAGS) -o libluna.so eval.o globals.o -Wl,--whole-archive *.a -Wl,--no-whole-archive
endif


libedf.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd edf; $(MAKE) $(MFLAGS)

libdb.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd db; $(MAKE) $(MFLAGS)

libsstore.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd sstore; $(MAKE) $(MFLAGS)

libgraphics.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd graphics; $(MAKE) $(MFLAGS)

libstaging.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd staging; $(MAKE) $(MFLAGS)

libdefs.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd defs; $(MAKE) $(MFLAGS)

libstats.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd stats; $(MAKE) $(MFLAGS)

libzfile.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd zfile; $(MAKE) $(MFLAGS)

libfftwrap.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd fftw; $(MAKE) $(MFLAGS)

libdsp.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd dsp; $(MAKE) $(MFLAGS)

libcwt.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd cwt; $(MAKE) $(MFLAGS)

libpdc.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd pdc; $(MAKE) $(MFLAGS)

libica.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd ica; $(MAKE) $(MFLAGS)

libclocs.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd clocs; $(MAKE) $(MFLAGS)

libartifacts.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd artifacts; $(MAKE) $(MFLAGS)

libspindles.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd spindles; $(MAKE) $(MFLAGS)

libtimeline.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd timeline; $(MAKE) $(MFLAGS)

libmiscmath.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd miscmath; $(MAKE) $(MFLAGS)

libannot.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd annot; $(MAKE) $(MFLAGS)

libintervals.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd intervals; $(MAKE) $(MFLAGS)

libhelper.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd helper; $(MAKE) $(MFLAGS)

libtinyxml.a : force_look
	$(ECHO) looking into subdir : $(MAKE) $(MFLAGS)
	cd tinyxml; $(MAKE) $(MFLAGS)

utils : force_look $(OBJLIBS)
	cd utils && $(MAKE)

clean :
	$(ECHO) cleaning up in .
	-$(RM) -f $(OBJS)
	-$(RM) -f *~
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done
	cd utils && $(MAKE) clean

cleanall :
	$(ECHO) cleaning up in .
	-$(RM) -f $(EXE) $(OBJS) $(OBJLIBS)
	-$(RM) -f *~
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done
	cd utils && $(MAKE) clean

force_look :
	true
