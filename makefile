CC=gcc
CFLAGS = -Wall -O2

all: GA

GA: GA.o functions.o OF.o
	$(CC) $(CFLAGS) -o GA GA.o functions.o OF.o -lm

GA.o: GA.c functions.h
	$(CC) $(CFLAGS) -c GA.c

functions.o: functions.c functions.h
	$(CC) $(CFLAGS) -c functions.c

OF.o: OF.c
	$(CC) $(CFLAGS) -c OF.c

clean:
	rm -f *.o GA