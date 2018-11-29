include ../Makefile.inc

CFLAGS += -Wpointer-sign

OBJLIBS	 = ../libclocs.a

OBJS = legendre_polynomial.o clocs.o topo.o

all : $(OBJLIBS)

../libclocs.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libclocs.a $?
	$(AR) $(ARFLAGS) rv ../libclocs.a $?
	$(ECHO) $(RANLIB) ../libclocs.a
	$(RANLIB) ../libclocs.a

clean :
	$(ECHO) cleaning up in clocs/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
