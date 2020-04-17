include ../Makefile.inc

OBJLIBS	 = ../libtimeline.a
OBJS	 = timeline.o hypno.o

all : $(OBJLIBS)

../libtimeline.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libtimeline.a $?
	$(AR) $(ARFLAGS) rv ../libtimeline.a $?
	$(ECHO) $(RANLIB) ../libtimeline.a
	$(RANLIB) ../libtimeline.a

clean :
	$(ECHO) cleaning up in timeline/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
