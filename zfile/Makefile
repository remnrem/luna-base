include ../Makefile.inc

OBJLIBS	 = ../libzfile.a
OBJS	 = zfstream.o 

all : $(OBJLIBS)

../libzfile.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libzfile.a $?
	$(AR) $(ARFLAGS) rv ../libzfile.a $?
	$(ECHO) $(RANLIB) ../libzfile.a
	$(RANLIB) ../libzfile.a

clean :
	$(ECHO) cleaning up in zfile/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
