#Sean M. Law

include ../makefile.defs

CFLAGS = -Wall $(LIBFLAGS)

LIBOBJS = cppmine.o mine.o 

all:	libmine.$(LIBEXT)

libmine.$(LIBEXT):	$(LIBOBJS)
	ar rc libmine.$(LIBEXT) $(LIBOBJS)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp  
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(DEBUG) -c $< -o $@
