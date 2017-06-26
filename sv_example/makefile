CXX       = g++
INITSRC   = ../base_model/initialize
ESTSRC    = ../base_model/estimator
ISCL      = ../base_model/libscl/gpp
LSCL      = ../base_model/libscl/gpp
USRSRC    = ./usermodel
CXXFLAGS  = -O2 -Wall -c -I$(INITSRC) -I$(USRSRC) -I$(ESTSRC) -I$(ISCL)
LDFLAGS   = -L$(LSCL) -lscl -lm

HEADERS  = $(ESTSRC)/estimator.h $(ESTSRC)/estimator_base.h $(INITSRC)/initialize.h \
	$(ISCL)/libscl.h $(INITSRC)/pathname.h \
	$(USRSRC)/usermodel.h $(USRSRC)/sv_moments.h $(USRSRC)/sv_model.h

OBJECTS = initialize.o mcmc_class.o proposal.o asymptotics.o \
	usermodel.o sv_moments.o sv_model.o

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

sv_moments.o : $(USRSRC)/sv_moments.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/sv_moments.cpp

sv_model.o : $(USRSRC)/sv_model.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/sv_model.cpp


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
