
CC  = gcc-mp-5
CXX = g++-mp-5
FC  = gfortran-mp-5

INSTALL_PATH?=/opt/local

SDIR = .
TRAME_DIR = $(SDIR)/../headers

CXX_STD = -std=c++11
OPT_FLAGS = -O2

SHLIB = libtrame.so
SHLIB_FLAGS = -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress

CXXFLAGS = -Wall $(OPT_FLAGS) -I/opt/local/include
FFLAGS = $(OPT_FLAGS)
LIBS= -lgfortran -framework Accelerate

LP_CFLAGS= -Wall $(OPT_FLAGS) -I/opt/local/include -I/Library/gurobi651/mac64/include
LP_LIBS= -L/Library/gurobi651/mac64/lib -lgurobi65
SOURCES_LP= $(TRAME_DIR)/lp/generic_lp_c.c
OBJECTS_LP= $(TRAME_DIR)/lp/generic_lp_c.o

LP_CXXFLAGS= -I/Library/gurobi651/mac64/include
SOURCES_LPXX= $(TRAME_DIR)/lp/generic_lp.cpp
OBJECTS_LPXX= $(TRAME_DIR)/lp/generic_lp.o

NLOPT_CXXFLAGS= -I/usr/local/include
NLOPT_LIBS = -L/usr/local/lib -lnlopt
SOURCES_NLOPT= $(TRAME_DIR)/nlopt/generic_nlopt.cpp
OBJECTS_NLOPT= $(TRAME_DIR)/nlopt/generic_nlopt.o

SOURCES_TRAME= $(TRAME_DIR)/arums/arums_empirical.cpp $(TRAME_DIR)/arums/arums_logit.cpp $(TRAME_DIR)/arums/arums_none.cpp $(TRAME_DIR)/arums/arums_probit.cpp $(TRAME_DIR)/arums/arums_rsc.cpp $(TRAME_DIR)/arums/arums_rusc.cpp
OBJECTS_TRAME= $(SOURCES_TRAME:.cpp=.o)

SOURCES= $(SDIR)/arums_logit_test.cpp $(SDIR)/arums_none_test.cpp $(SDIR)/arums_probit_test.cpp $(SDIR)/arums_rsc_test.cpp
OBJECTS= $(SOURCES:.cpp=.o)
PROGS= $(SOURCES:.cpp=.test)

OBJECTS_AUX= $(TRAME_DIR)/aux/inv_pwa.o $(TRAME_DIR)/aux/trame_aux.o $(TRAME_DIR)/aux/trame_stats.o $(TRAME_DIR)/aux/zeroin.o $(TRAME_DIR)/prob/prob.o $(TRAME_DIR)/math/quadpack_double.o $(TRAME_DIR)/prob/trame_F90_aux.o

all: $(PROGS) $(TRAME_DIR)/$(SHLIB) $(OBJECTS) $(OBJECTS_TRAME) $(OBJECTS_LPXX) $(OBJECTS_NLOPT) $(OBJECTS_AUX)

#$(SDIR)/%.test: $(SDIR)/%.o | $(OBJECTS_LP) $(OBJECTS_LPXX) $(OBJECTS_NLOPT) $(OBJECTS_AUX) $(OBJECTS_TRAME)
#	$(CXX) -Wall $(OPT_FLAGS) $(CXX_STD) -o $@ $< $(OBJECTS_LP) $(OBJECTS_LPXX) $(OBJECTS_NLOPT) $(OBJECTS_AUX) $(OBJECTS_TRAME) $(LP_LIBS) $(NLOPT_LIBS) $(LIBS)

$(SDIR)/%.test: $(SDIR)/%.o | $(TRAME_DIR)/$(SHLIB)
	$(CXX) -Wall $(OPT_FLAGS) $(CXX_STD) -o $@ $< -L$(TRAME_DIR) -ltrame -framework Accelerate

$(TRAME_DIR)/$(SHLIB): $(OBJECTS_LP) $(OBJECTS_LPXX) $(OBJECTS_NLOPT) $(OBJECTS_AUX) $(OBJECTS_TRAME)
	$(CXX) $(CXX_STD) $(SHLIB_FLAGS) -o $@ $^ $(LP_LIBS) $(NLOPT_LIBS) $(LIBS)

$(SDIR)/%.o: $(SDIR)/%.cpp
	$(CXX) $(CXX_STD) $(CXXFLAGS) $(LP_CXXFLAGS) $(NLOPT_CXXFLAGS) $< -c -o $@ 

#$(OBJECTS): $(SOURCES) $(OBJECTS_LP) $(OBJECTS_LPXX) $(OBJECTS_NLOPT) $(OBJECTS_AUX) $(OBJECTS_TRAME)
#	$(CXX) $(CXX_STD) $(CXXFLAGS) $< -o $< $(LP_LIBS) $(NLOPT_LIBS) $(LIBS)

$(TRAME_DIR)/arums/%.o: $(TRAME_DIR)/arums/%.cpp
	$(CXX) $(CXX_STD) $(CXXFLAGS) $(LP_CXXFLAGS) $(NLOPT_CXXFLAGS) $< -c -o $@

$(OBJECTS_LP): $(SOURCES_LP)
	$(CC) $(LP_CFLAGS) $^ -c -o $@

$(OBJECTS_LPXX): $(SOURCES_LPXX)
	$(CXX) $(CXXFLAGS) $(LP_CXXFLAGS) $^ -c -o $@

$(OBJECTS_NLOPT): $(SOURCES_NLOPT)
	$(CXX) $(CXXFLAGS) $(NLOPT_CXXFLAGS) $^ -c -o $@

$(TRAME_DIR)/aux/%.o: $(TRAME_DIR)/aux/%.cpp
	$(CXX) $(CXX_STD) $(CXXFLAGS) $< -c -o $@

#$(TRAME_DIR)/aux/trame_aux.o: $(TRAME_DIR)/aux/trame_aux.cpp
#	$(CXX) $(CXXFLAGS) $^ -c -o $@

$(TRAME_DIR)/prob/prob.o: $(TRAME_DIR)/prob/prob.f90
	$(FC) $(FFLAGS) $^ -c -o $@

$(TRAME_DIR)/math/quadpack_double.o: $(TRAME_DIR)/math/quadpack_double.f90
	$(FC) $(FFLAGS) $^ -c -o $@

$(TRAME_DIR)/prob/trame_F90_aux.o: $(TRAME_DIR)/prob/trame_F90_aux.f90
	$(FC) $(FFLAGS) $^ -c -o $@

clean:
	rm -f ./*.o ./*.test $(TRAME_DIR)/arums/*.o $(TRAME_DIR)/aux/*.o $(TRAME_DIR)/lp/*.o $(TRAME_DIR)/math/*.o $(TRAME_DIR)/nlopt/*.o $(TRAME_DIR)/prob/*.o

.PHONY: install
install: $(SHLIB)
	#mkdir $(INSTALL_PATH)/include
	#mkdir $(INSTALL_PATH)/lib
	cp $(TRAME_DIR)/$(SHLIB) $(INSTALL_PATH)/lib/$(SHLIB)
	cp $(TRAME_DIR)/trame.hpp $(INSTALL_PATH)/include

