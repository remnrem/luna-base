include ../Makefile.inc

OBJLIBS	 = ../libhelper.a
OBJS	 = helper.o

all : $(OBJLIBS)

../libhelper.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libhelper.a $?
	$(AR) $(ARFLAGS) rv ../libhelper.a $?
	$(ECHO) $(RANLIB) ../libhelper.a
	$(RANLIB) ../libhelper.a

clean :
	$(ECHO) cleaning up in helper/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
