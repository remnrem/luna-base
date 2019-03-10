include ../../Makefile.inc

CXXFLAGS += -I../../

OBJLIB	 = ../../libsrate.a

OBJS = samplerate.o src_linear.o src_sinc.o src_zoh.o

all : $(OBJLIB)

../../libsrate.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv $(OBJLIB) $?
	$(AR) $(ARFLAGS) rv $(OBJLIB) $?
	$(ECHO) $(RANLIB) $(OBJLIB)
	$(RANLIB) $(OBJLIB)

clean :
	$(ECHO) cleaning up in libsamplerate/
	-$(RM) -f $(OBJS) 
	-$(RM) -f *~

force_look :
	true
