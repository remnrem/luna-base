include ../Makefile.inc

OBJLIBS	 = ../libartifacts.a
OBJS	 = artifacts.o

all : $(OBJLIBS)

../libartifacts.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libartifacts.a $?
	$(AR) $(ARFLAGS) rv ../libartifacts.a $?
	$(ECHO) $(RANLIB) ../libartifacts.a
	$(RANLIB) ../libartifacts.a

clean :
	$(ECHO) cleaning up in artifacts/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
