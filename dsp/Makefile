include ../Makefile.inc

CFLAGS += -Wpointer-sign

OBJLIBS	 = ../libdsp.a

OBJS = reduce.o mi.o resample.o coherence.o pac.o phsyn.o ecgsuppression.o hilbert.o	\
slow-waves.o emd.o mse.o cfc.o lzw.o  fir.o fiplot.o ed.o	\
interpolate.o correl.o conv.o polarity.o spectral_norm.o cwt-design.o 	\
r8lib.o pwl_interp_2d_scattered.o tv.o wrappers.o ica-wrapper.o sl.o

all : $(OBJLIBS)

../libdsp.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libdsp.a $?
	$(AR) $(ARFLAGS) rv ../libdsp.a $?
	$(ECHO) $(RANLIB) ../libdsp.a
	$(RANLIB) ../libdsp.a

clean :
	$(ECHO) cleaning up in dsp/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
