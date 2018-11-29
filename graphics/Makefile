include ../Makefile.inc

OBJLIBS	 = ../libgraphics.a
OBJS	 = graphics.o

all : $(OBJLIBS)

../libgraphics.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libgraphics.a $?
	$(AR) $(ARFLAGS) rv ../libgraphics.a $?
	$(ECHO) $(RANLIB) ../libgraphics.a
	$(RANLIB) ../libgraphics.a

clean :
	$(ECHO) cleaning up in graphics/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
