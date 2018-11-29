include ../Makefile.inc

OBJLIBS	 = ../libspindles.a
OBJS	 = spindles.o mspindles.o plot-spindles.o spectral.o

all : $(OBJLIBS)

../libspindles.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libspindles.a $?
	$(AR) $(ARFLAGS) rv ../libspindles.a $?
	$(ECHO) $(RANLIB) ../libspindles.a
	$(RANLIB) ../libspindles.a

clean :
	$(ECHO) cleaning up in spindles/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
