include ../Makefile.inc

OBJLIBS	 = ../libsuds.a
OBJS	 = suds.o

all : $(OBJLIBS)

../libsuds.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libsuds.a $?
	$(AR) $(ARFLAGS) rv ../libsuds.a $?
	$(ECHO) $(RANLIB) ../libsuds.a
	$(RANLIB) ../libsuds.a

clean :
	$(ECHO) cleaning up in suds/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
