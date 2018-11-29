include ../Makefile.inc

OBJLIBS	 = ../libcwt.a
OBJS	 = cwt.o

all : $(OBJLIBS)

../libcwt.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libcwt.a $?
	$(AR) $(ARFLAGS) rv ../libcwt.a $?
	$(ECHO) $(RANLIB) ../libcwt.a
	$(RANLIB) ../libcwt.a

clean :
	$(ECHO) cleaning up in cwt/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
