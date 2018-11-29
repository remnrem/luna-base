include ../Makefile.inc

OBJLIBS	 = ../libintervals.a
OBJS	 = intervals.o

all : $(OBJLIBS)

../libintervals.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libintervals.a $?
	$(AR) $(ARFLAGS) rv ../libintervals.a $?
	$(ECHO) $(RANLIB) ../libintervals.a
	$(RANLIB) ../libintervals.a

clean :
	$(ECHO) cleaning up in intervals/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
