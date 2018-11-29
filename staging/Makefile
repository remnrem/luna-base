include ../Makefile.inc

OBJLIBS	 = ../libstaging.a
OBJS	 = zratio.o

all : $(OBJLIBS)

../libstaging.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libstaging.a $?
	$(AR) $(ARFLAGS) rv ../libstaging.a $?
	$(ECHO) $(RANLIB) ../libstaging.a
	$(RANLIB) ../libstaging.a

clean :
	$(ECHO) cleaning up in staging/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
