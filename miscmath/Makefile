include ../Makefile.inc

OBJLIBS	 = ../libmiscmath.a
OBJS	 = miscmath.o dynam.o crandom.o

all : $(OBJLIBS)

../libmiscmath.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libmiscmath.a $?
	$(AR) $(ARFLAGS) rv ../libmiscmath.a $?
	$(ECHO) $(RANLIB) ../libmiscmath.a
	$(RANLIB) ../libmiscmath.a

clean :
	$(ECHO) cleaning up in miscmath/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
