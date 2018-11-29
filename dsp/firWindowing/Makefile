CC			= g++
CFLAGS	+= -Wall -O2
LIBRARY	= -lm -lfftw3

all: filters

filters: filters.c
	$(CC) $(CFLAGS) $(LIBRARY) -o filters filters.cpp

