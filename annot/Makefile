include ../Makefile.inc

OBJLIBS	 = ../libannot.a
OBJS	 = annot.o

all : $(OBJLIBS)

../libannot.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libannot.a $?
	$(AR) $(ARFLAGS) rv ../libannot.a $?
	$(ECHO) $(RANLIB) ../libannot.a
	$(RANLIB) ../libannot.a

clean :
	$(ECHO) cleaning up in annot/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
