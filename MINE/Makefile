#Sean M. Law

include ../makefile.defs

CFLAGS = -Wall -lm $(LIBFLAGS)

LIBOBJS = cppmine.o mine.o 

all:	libmine.$(LIBEXT)

libmine.$(LIBEXT):	$(LIBOBJS)
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $^ -o libmine.$(LIBEXT)

clean:
	rm -rf *.o
	rm -rf *.$(LIBEXT)
	rm -rf *stackdump
	rm -rf *~

%.o:	%.cpp  
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) $(LIBFLAGS) $(DEBUG) -c $< -o $@
