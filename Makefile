#

CC=g++
CFLAGS=-Wall -Ilib

all:
	mkdir -p bin
	$(CC) $(CFLAGS) src/manipPDB.c -o bin/manipPDB 

clean:
	rm -rf bin *~ src/*~ lib/*~ tests/*~
