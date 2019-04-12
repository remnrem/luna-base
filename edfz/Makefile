include ../Makefile.inc

OBJLIBS	 = ../libedfz.a
OBJS	 = bgzf.o edfz.o

all : $(OBJLIBS)

../libedfz.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libedfz.a $?
	$(AR) $(ARFLAGS) rv ../libedfz.a $?
	$(ECHO) $(RANLIB) ../libedfz.a
	$(RANLIB) ../libedfz.a

clean :
	$(ECHO) cleaning up in edfz/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
