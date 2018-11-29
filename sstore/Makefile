include ../Makefile.inc

OBJLIBS	 = ../libsstore.a
OBJS	 = ../db/sqlite3.o ../db/sqlwrap.o sstore.o

all : $(OBJLIBS)

../libsstore.a : $(OBJS)
	$(ECHO) $(AR) $(ARFLAGS) rv ../libsstore.a $?
	$(AR) $(ARFLAGS) rv ../libsstore.a $?
	$(ECHO) $(RANLIB) ../libsstore.a
	$(RANLIB) ../libsstore.a

clean :
	$(ECHO) cleaning up in sstore/
	-$(RM) -f $(OBJS)
	-$(RM) -f *~

force_look :
	true
