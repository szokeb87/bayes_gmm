CXX       = g++
ISCL      = ../libscl/gpp
LSCL      = ../libscl/gpp
CXXFLAGS  = -O2 -Wall -c -I$(ISCL)
LDFLAGS   = -L$(LSCL) -lscl -lm

HEADERS  = $(ISCL)/libscl.h

LIBRARIES = libscl.a

PROGRAMS = bayes_gmm

all : $(LIBRARIES) $(PROGRAMS)

libscl.a :
	make -C $(LSCL)

bayes_gmm :  test.o $(HEADERS) $(LIBRARIES)
	$(CXX) -o bayes_gmm test.o $(LDFLAGS)
	rm -f *.o core core.*

test.o : ./test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) ./test.cpp


clean	:
	make -C $(LSCL) clean
	rm -f *.o core core.*

veryclean	:
	make -C $(LSCL) veryclean
	rm -f *.o core core.*
	rm -f $(PROGRAMS)

dirclean	:
	rm -f *.o core core.*
	rm -f $(PROGRAMS)
