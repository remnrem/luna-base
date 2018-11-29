include ../Makefile.inc

OBJLIBS	 = ../libdb.a
OBJS	 = sqlite3.o db.o sqlwrap.o

all : $(OBJLIBS)

../libdb.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libdb.a $?
	$(AR) $(ARFLAGS) rv ../libdb.a $?
	$(ECHO) $(RANLIB) ../libdb.a
	$(RANLIB) ../libdb.a

clean :
	$(ECHO) cleaning up in db/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
