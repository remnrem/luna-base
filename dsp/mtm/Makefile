include ../../Makefile.inc

CXXFLAGS += -I../../

OBJLIB	 = ../../libmtm.a

OBJS = multitap.o jtridib.o jtinvit.o jrealft.o jfour1.o hires.o	\
ftest.o dfour1.o adwait.o sigstuff.o nrutil.o mult_tap_spec.o mtm.o

all : $(OBJLIB)

../../libmtm.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv $(OBJLIB) $?
	$(AR) $(ARFLAGS) rv $(OBJLIB) $?
	$(ECHO) $(RANLIB) $(OBJLIB)
	$(RANLIB) $(OBJLIB)

clean :
	$(ECHO) cleaning up in mtm/
	-$(RM) -f $(OBJS) 
	-$(RM) -f *~

force_look :
	true
