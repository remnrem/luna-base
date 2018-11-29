include ../Makefile.inc

CFLAGS += -Wpointer-sign

OBJLIBS	 = ../libica.a

OBJS = ica.o libICA.o matrix.o svdcmp.o

all : $(OBJLIBS)

../libica.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libica.a $?
	$(AR) $(ARFLAGS) rv ../libica.a $?
	$(ECHO) $(RANLIB) ../libica.a
	$(RANLIB) ../libica.a

clean :
	$(ECHO) cleaning up in ica/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
