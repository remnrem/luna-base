include ../Makefile.inc

OBJLIBS	 = ../libdefs.a
OBJS	 = defs.o

all : $(OBJLIBS)

../libdefs.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libdefs.a $?
	$(AR) $(ARFLAGS) rv ../libdefs.a $?
	$(ECHO) $(RANLIB) ../libdefs.a
	$(RANLIB) ../libdefs.a

clean :
	$(ECHO) cleaning up in defs/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
