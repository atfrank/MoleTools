#

CC=g++
CFLAGS=-Wall -Ilib

all:
	mkdir -p bin
	$(CC) $(CFLAGS) src/manipPDB.cpp -o bin/manipPDB 

clean:
	rm -rf bin *~ src/*~ lib/*~ tests/*~
