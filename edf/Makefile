include ../Makefile.inc

OBJLIBS	 = ../libedf.a
OBJS	 = edf.o slice.o dumper.o covar.o tal.o masks.o dump-intervals.o chep.o
 
all : $(OBJLIBS)

../libedf.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libedf.a $?
	$(AR) $(ARFLAGS) rv ../libedf.a $?
	$(ECHO) $(RANLIB) ../libedf.a
	$(RANLIB) ../libedf.a

clean :
	$(ECHO) cleaning up in edf/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
