include Makefile.inc

DIRS = edf tinyxml helper timeline annot dsp miscmath spindles	\
artifacts intervals fftw cwt defs stats graphics staging 	\
db ica clocs pdc sstore dsp/mtm dsp/libsamplerate

EXE	= luna
LUNALIB = libluna.so
OBJS	= main.o globals.o eval.o

OBJLIBS = libdefs.a libedf.a libtinyxml.a libhelper.a libtimeline.a	\
libannot.a libdsp.a libmiscmath.a libspindles.a libartifacts.a		\
libintervals.a libfftwrap.a libcwt.a libstats.a libgraphics.a		\
libstaging.a libdb.a libica.a libclocs.a libpdc.a libsstore.a libmtm.a	\
libsrate.a

LIBS = -L. -lspindles -lica -lannot -ldefs -lartifacts -ledf -lhelper	\
-ltimeline -lstaging -lfftwrap -ldsp -lmtm -lmiscmath -lintervals	\
-ltinyxml -lcwt -lclocs -lpdc -lstats -lgraphics -ldb -lsstore -lsrate	\
-lfftw3

all : $(EXE) $(LUNALIB) utils

$(EXE) : main.o globals.o eval.o $(OBJLIBS)
	$(ECHO) $(LD) $(LDFLAGS) -o $(EXE) $(OBJS) $(LIBS)
	$(LD) $(LDFLAGS) -o $(EXE) $(OBJS) $(LIBS)

static : main.o globals.o eval.o $(OBJLIBS)
	g++ -static -static-libgcc -static-libstdc++ -L/usr/local/lib	\
	-o $(LUNA_STATIC) main.o globals.o eval.o libspindles.a libartifacts.a	\
	libtimeline.a libannot.a libedf.a libintervals.a libcwt.a libdsp.a	\
	libstaging.a libclocs.a libpdc.a libsstore.a libmtm.a libdefs.a		\
	libhelper.a libfftwrap.a libgraphics.a libmiscmath.a libstats.a		\
	libsrate.a libtinyxml.a libdb.a libica.a $(FFTW)/lib/libfftw3.a


$(LUNALIB) : globals.o eval.o $(OBJLIBS)
ifeq ($(ARCH),MAC)
	$(ECHO) "building libluna.dylib..."
	$(LD) -dynamiclib -fPIC $(LDFLAGS) -o libluna.dylib eval.o globals.o  *.a  -lfftw3
else
	$(ECHO) "building libluna.so..."
	$(LD) -shared     -fPIC $(LDFLAGS) -o libluna.so eval.o globals.o -Wl,--whole-archive *.a -Wl,--no-whole-archive
endif


libedf.a : force_look
	cd edf; $(MAKE) $(MFLAGS)

libdb.a : force_look
	cd db; $(MAKE) $(MFLAGS)

libsstore.a : force_look
	cd sstore; $(MAKE) $(MFLAGS)

libgraphics.a : force_look
	cd graphics; $(MAKE) $(MFLAGS)

libstaging.a : force_look
	cd staging; $(MAKE) $(MFLAGS)

libdefs.a : force_look
	cd defs; $(MAKE) $(MFLAGS)

libstats.a : force_look
	cd stats; $(MAKE) $(MFLAGS)

#libzfile.a : force_look
#	cd zfile; $(MAKE) $(MFLAGS)

libfftwrap.a : force_look
	cd fftw; $(MAKE) $(MFLAGS)

libdsp.a : force_look
	cd dsp; $(MAKE) $(MFLAGS)

libmtm.a : force_look
	cd dsp/mtm; $(MAKE) $(MFLAGS)

libsrate.a : force_look
	cd dsp/libsamplerate; $(MAKE) $(MFLAGS)

libcwt.a : force_look
	cd cwt; $(MAKE) $(MFLAGS)

libpdc.a : force_look
	cd pdc; $(MAKE) $(MFLAGS)

libica.a : force_look
	cd ica; $(MAKE) $(MFLAGS)

libclocs.a : force_look
	cd clocs; $(MAKE) $(MFLAGS)

libartifacts.a : force_look
	cd artifacts; $(MAKE) $(MFLAGS)

libspindles.a : force_look
	cd spindles; $(MAKE) $(MFLAGS)

libtimeline.a : force_look
	cd timeline; $(MAKE) $(MFLAGS)

libmiscmath.a : force_look
	cd miscmath; $(MAKE) $(MFLAGS)

libannot.a : force_look
	cd annot; $(MAKE) $(MFLAGS)

libintervals.a : force_look
	cd intervals; $(MAKE) $(MFLAGS)

libhelper.a : force_look
	cd helper; $(MAKE) $(MFLAGS)

libtinyxml.a : force_look
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
