include ./makefile.defs

all: libmoletools.$(LIBEXT) libmine.$(LIBEXT) programs

libmoletools.$(LIBEXT):
	$(UNAME)
	$(MAKE) -C MoleTools 

libmine.$(LIBEXT):
	$(UNAME)
	$(MAKE) -C MINE

programs:
	$(MAKE) -C src

clean:
	$(MAKE) -C MoleTools clean
	$(MAKE) -C MINE clean
	$(MAKE) -C src clean
	rm -rf *~
	rm -rf *stackdump
