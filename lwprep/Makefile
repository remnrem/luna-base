include ../Makefile.inc

OBJLIBS	 = ../liblwprep.a
OBJS	 = lwprep.o

all : $(OBJLIBS)

../liblwprep.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../liblwprep.a $?
	$(AR) $(ARFLAGS) rv ../liblwprep.a $?
	$(ECHO) $(RANLIB) ../liblwprep.a
	$(RANLIB) ../liblwprep.a

clean :
	$(ECHO) cleaning up in lwprep/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
