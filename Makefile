#

CC=g++
CFLAGS=-c -Wall -Ilib

all:
	mkdir bin
	$(CC) $(CFLAGS) src/manipPDB.c -o bin/manipPDB

clean:
	rm -rf bin *~ src/*~ lib/*~
