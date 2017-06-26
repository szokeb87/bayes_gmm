CXX       = g++
USRSRC    = ../usermodel
ISCL      = ../../base_model/libscl/gpp
LSCL      = ../../base_model/libscl/gpp
CXXFLAGS  = -O2 -Wall -c -I$(USRSRC) -I$(ISCL)
LDFLAGS   = -L$(LSCL) -lscl -lm

HEADERS  = $(ISCL)/libscl.h $(USRSRC)/sv_model.h $(USRSRC)/sv_moments.h

LIBRARIES = libscl.a

PROGRAMS = sv_simul_data

all :  $(LIBRARIES) $(PROGRAMS)

libscl.a :
	make -C $(LSCL)

sv_simul_data :  sv_simul_data.o sv_model.o	sv_moments.o $(HEADERS) $(LIBRARIES)
	$(CXX) -o sv_simul_data sv_simul_data.o sv_model.o sv_moments.o $(LDFLAGS)
	rm -f *.o core core.*

sv_simul_data.o : sv_simul_data.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) sv_simul_data.cpp

sv_model.o : $(USRSRC)/sv_model.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/sv_model.cpp

sv_moments.o : $(USRSRC)/sv_moments.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(USRSRC)/sv_moments.cpp
