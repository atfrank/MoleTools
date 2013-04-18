include ./makefile.defs

all: libtools.so programs

libtools.so:
	$(MAKE) -C lib 

programs:
	$(MAKE) -C src

clean:
	$(MAKE) -C lib clean
	$(MAKE) -C src clean
