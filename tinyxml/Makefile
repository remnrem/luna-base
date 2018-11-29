include ../Makefile.inc

OBJLIBS	 = ../libtinyxml.a
OBJS	 = tinyxml.o tinyxmlerror.o tinyxmlparser.o tinystr.o xmlreader.o

all : $(OBJLIBS)

../libtinyxml.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libtinyxml.a $?
	$(AR) $(ARFLAGS) rv ../libtinyxml.a $?
	$(ECHO) $(RANLIB) ../libtinyxml.a
	$(RANLIB) ../libtinyxml.a

clean :
	$(ECHO) cleaning up in tinyxml/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
