include ../Makefile.inc

OBJLIBS	 = ../libstats.a
OBJS	 = dcdflib.o matrix.o statistics.o glm.o cluster.o


all : $(OBJLIBS)

../libstats.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libstats.a $?
	$(AR) $(ARFLAGS) rv ../libstats.a $?
	$(ECHO) $(RANLIB) ../libstats.a
	$(RANLIB) ../libstats.a

clean :
	$(ECHO) cleaning up in stats/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
