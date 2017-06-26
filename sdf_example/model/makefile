CXX       = g++
INITSRC   = ../model/initialize
ESTSRC    = ../model/estimator
USRSRC    = ../model/usermodel
ISCL      = ../model/libscl/gpp
LSCL      = ../model/libscl/gpp
CXXFLAGS  = -O2 -Wall -c -I$(INITSRC) -I$(USRSRC) -I$(ESTSRC) -I$(ISCL)
LDFLAGS   = -L$(LSCL) -lscl -lm

HEADERS  = $(ESTSRC)/estimator.h $(ESTSRC)/estimator_base.h \
	$(INITSRC)/initialize.h $(INITSRC)/pathname.h \
	$(ISCL)/libscl.h \
	$(USRSRC)/usermodel.h $(USRSRC)/moments.h $(USRSRC)/model.h $(USRSRC)/default_params.h

OBJECTS = initialize.o mcmc_class.o proposal.o asymptotics.o \
	usermodel.o moments.o model.o default_params.o

LIBRARIES = libscl.a

PROGRAMS = bayes_gmm

all : pathname.h $(LIBRARIES) $(PROGRAMS)

libscl.a :
	make -C $(LSCL)

bayes_gmm :  main.o $(OBJECTS) $(HEADERS) $(LIBRARIES)
	$(CXX) -o bayes_gmm main.o $(OBJECTS) $(LDFLAGS)
	rm -f *.o core core.*

pathname.h:
	echo "#define PATHNAME \"$(PWD)\"" > $(INITSRC)/pathname.h

main.o : $(INITSRC)/main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INITSRC)/main.cpp

initialize.o : $(INITSRC)/initialize.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INITSRC)/initialize.cpp

mcmc_class.o : $(ESTSRC)/mcmc_class.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ESTSRC)/mcmc_class.cpp

proposal.o : $(ESTSRC)/proposal.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ESTSRC)/proposal.cpp

asymptotics.o : $(ESTSRC)/asymptotics.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(ESTSRC)/asymptotics.cpp

usermodel.o : $(USRSRC)/usermodel.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/usermodel.cpp

moments.o : $(USRSRC)/moments.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/moments.cpp

model.o : $(USRSRC)/model.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/model.cpp

default_params.o : $(USRSRC)/default_params.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/default_params.cpp

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
