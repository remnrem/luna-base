
make
g++ -O2 -c loadss.cpp
g++ -O2 -o loadss loadss.o sstore.o ../db/sqlwrap.o ../db/sqlite3.o ../helper/helper.o ../defs/defs.o ../miscmath/crandom.o

# simple interval loader :  loadints name < {ints}  
g++ -O2 -c loadints.cpp
g++ -O2 -o loadints loadints.o sstore.o ../db/sqlwrap.o ../db/sqlite3.o ../helper/helper.o ../defs/defs.o ../miscmath/crandom.o

g++ -O2 -c tabless.cpp
g++ -O2 -o tabless tabless.o sstore.o ../db/sqlwrap.o ../db/sqlite3.o ../helper/helper.o ../defs/defs.o ../miscmath/crandom.o
