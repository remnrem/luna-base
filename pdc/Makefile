include ../Makefile.inc

OBJLIBS	 = ../libpdc.a
OBJS	 = pdc.o pdcfuncs.o exe.o sss.o external.o

all : $(OBJLIBS)

../libpdc.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libpdc.a $?
	$(AR) $(ARFLAGS) rv ../libpdc.a $?
	$(ECHO) $(RANLIB) ../libpdc.a
	$(RANLIB) ../libpdc.a

clean :
	$(ECHO) cleaning up in pdc/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
