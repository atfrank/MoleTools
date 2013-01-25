#

CC=g++
CFLAGS=-c -Wall

all:
	mkdir bin
	$(CC) $(CFLAGS) src/manipPDB.c -o bin/manipPDB

clean:
	rm -rf bin
