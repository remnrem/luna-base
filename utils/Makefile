include ../Makefile.inc

# set all executables to be placed in parent directory

DESTRAT = ../destrat
DESTRAT_OBJS = reader.o ../globals.o
DESTRAT_LIBS = -L.. -lhelper -ldefs -lmiscmath -ldb -lannot

TSPICKER = ../tspicker 
TSPICKER_OBJS = tspicker.o ../globals.o
TSPICKER_LIBS = -L.. -lhelper -ldefs -lmiscmath -ldb -lannot

INTERSECT = ../intersect 
INTERSECT_OBJS = list-intersection.o ../globals.o
INTERSECT_LIBS = -L.. -lhelper -ldefs -lmiscmath -lintervals -ldb -lannot

all : $(DESTRAT) $(TSPICKER) $(INTERSECT)

$(DESTRAT) : reader.o ../globals.o
	$(ECHO) $(LD) $(LDFLAGS) -o $(DESTRAT) $(DESTRAT_OBJS) $(DESTRAT_LIBS)
	$(LD) $(LDFLAGS) -o $(DESTRAT) $(DESTRAT_OBJS) $(DESTRAT_LIBS)


$(TSPICKER) : tspicker.o ../globals.o
	$(ECHO) $(LD) $(LDFLAGS) -o $(TSPICKER) $(TSPICKER_OBJS) $(TSPICKER_LIBS)
	$(LD) $(LDFLAGS) -o $(TSPICKER) $(TSPICKER_OBJS) $(TSPICKER_LIBS)

$(INTERSECT) : list-intersection.o ../globals.o
	$(ECHO) $(LD) $(LDFLAGS) -o $(INTERSECT) $(INTERSECT_OBJS) $(INTERSECT_LIBS)
	$(LD) $(LDFLAGS) -o $(INTERSECT) $(INTERSECT_OBJS) $(INTERSECT_LIBS)

clean :
	$(ECHO) cleaning up utils/
	-$(RM) -f $(DESTRAT) $(DESTRAT_OBJS)
	-$(RM) -f $(TSPICKER) $(TSPICKER_OBJS) 
	-$(RM) -f $(INTERSECT) $(INTERSECT_OBJS) 
	-$(RM) -f *~
