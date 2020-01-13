include ../Makefile.inc

OBJLIBS	 = ../libfftwrap.a
OBJS	 = fftwrap.o cohfft.o

all : $(OBJLIBS)

../libfftwrap.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libfftwrap.a $?
	$(AR) $(ARFLAGS) rv ../libfftwrap.a $?
	$(ECHO) $(RANLIB) ../libfftwrap.a
	$(RANLIB) ../libfftwrap.a

clean :
	$(ECHO) cleaning up in fftwrap/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
