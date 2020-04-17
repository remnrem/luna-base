include ../Makefile.inc

# set all executables to be placed in parent directory

DESTRAT = ../destrat
DESTRAT_OBJS = reader.o ../globals.o ../cmddefs.o
DESTRAT_LIBS = -L..  -ldb -lsstore -lhelper -ldefs -lannot -lmiscmath -lz

BEHEAD = ../behead
BEHEAD_OBJS = behead.o
BEHEAD_LIBS = 

TSPICKER = ../tspicker 
TSPICKER_OBJS = tspicker.o ../globals.o ../cmddefs.o
TSPICKER_LIBS = -L.. -ldb -lsstore -lhelper -ldefs -lmiscmath -lannot  -lz

INTERSECT = ../intersect 
INTERSECT_OBJS = list-intersection.o ../globals.o ../cmddefs.o
INTERSECT_LIBS = -L.. -ldb -lsstore  -lhelper -ldefs -lmiscmath -lintervals -lannot  -lz

PREPSS = ../prepss
PREPSS_OBJS = prepss.o ../globals.o ../cmddefs.o
PREPSS_LIBS = -L.. -ldb -lsstore  -lhelper -ldefs -lmiscmath -lintervals -lannot  -lz

LOADSS = ../loadss
LOADSS_OBJS = loadss.o ../globals.o ../cmddefs.o
LOADSS_LIBS = -L.. -ldb -lsstore  -lhelper -ldefs -lmiscmath -lintervals -lannot  -lz

MERGE = ../merge
MERGE_OBJS = merge.o merge-helpers.o
MERGE_LIBS = 


all : $(DESTRAT) $(TSPICKER) $(INTERSECT) $(BEHEAD) $(PREPSS) $(LOADSS) $(MERGE)

$(DESTRAT) : ${DESTRAT_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(DESTRAT) $(DESTRAT_OBJS) $(DESTRAT_LIBS)
	$(LD) $(LDFLAGS) -o $(DESTRAT) $(DESTRAT_OBJS) $(DESTRAT_LIBS)

$(MERGE) : ${MERGE_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(MERGE) $(MERGE_OBJS) $(MERGE_LIBS)
	$(LD) $(LDFLAGS) -o $(MERGE) $(MERGE_OBJS) $(MERGE_LIBS)

$(BEHEAD) : ${BEHEAD_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(BEHEAD) $(BEHEAD_OBJS) $(BEHEAD_LIBS)
	$(LD) $(LDFLAGS) -o $(BEHEAD) $(BEHEAD_OBJS) $(BEHEAD_LIBS)

$(TSPICKER) : ${TSPICKER_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(TSPICKER) $(TSPICKER_OBJS) $(TSPICKER_LIBS)
	$(LD) $(LDFLAGS) -o $(TSPICKER) $(TSPICKER_OBJS) $(TSPICKER_LIBS)

$(INTERSECT) : $(INTERSECT_OBJS)
	$(ECHO) $(LD) $(LDFLAGS) -o $(INTERSECT) $(INTERSECT_OBJS) $(INTERSECT_LIBS)
	$(LD) $(LDFLAGS) -o $(INTERSECT) $(INTERSECT_OBJS) $(INTERSECT_LIBS)

$(PREPSS) : ${PREPSS_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(PREPSS) $(PREPSS_OBJS) $(PREPSS_LIBS)
	$(LD) $(LDFLAGS) -o $(PREPSS) $(PREPSS_OBJS) $(PREPSS_LIBS)

$(LOADSS) : ${LOADSS_OBJS}
	$(ECHO) $(LD) $(LDFLAGS) -o $(LOADSS) $(LOADSS_OBJS) $(LOADSS_LIBS)
	$(LD) $(LDFLAGS) -o $(LOADSS) $(LOADSS_OBJS) $(LOADSS_LIBS)

clean :
	$(ECHO) cleaning up utils/
	-$(RM) -f $(DESTRAT) $(DESTRAT_OBJS)
	-$(RM) -f $(BEHEAD) $(BEHEAD_OBJS)
	-$(RM) -f $(TSPICKER) $(TSPICKER_OBJS)
	-$(RM) -f $(INTERSECT) $(INTERSECT_OBJS)
	-$(RM) -f $(PREPSS) $(PREPSS_OBJS)
	-$(RM) -f $(LOADSS) $(LOADSS_OBJS)
	-$(RM) -f $(MERGE) $(MERGE_OBJS)
	-$(RM) -f *~
