
# core compiling options
#CC  = gcc-mp-5
#CXX = g++-mp-5
#FC  = gfortran-mp-5
CC  = gcc-mp-5
CXX = g++-mp-5
FC  = gfortran-mp-5

CXX_STD = -std=c++11
OPT_FLAGS = -O2

ARMA_INCLUDE_PATH = -I/opt/local/include

# install location
INSTALL_PATH=/opt/local

# source directories
SDIR = .
TRAME_DIR = $(SDIR)
TRAME_SRC_DIR = $(SDIR)/src
TRAME_HEADER_DIR = $(SDIR)/include
TRAME_TEST_DIR = $(SDIR)/tests

# shared library name and flags
SHLIB = libtrame.so
SHLIB_FLAGS = $(CXX_STD) -dynamiclib -install_name /opt/local/lib/libtrame.so -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress

# linear programming flags
LP_CFLAGS= -DTRAME_USE_GUROBI -I/Library/gurobi701/mac64/include
LP_LIBS= -L/Library/gurobi701/mac64/lib -lgurobi70
SOURCES_LP= $(TRAME_SRC_DIR)/lp/generic_lp_c.c
OBJECTS_LP= $(TRAME_SRC_DIR)/lp/generic_lp_c.o

LP_CXXFLAGS= -DTRAME_USE_GUROBI -I/Library/gurobi701/mac64/include
SOURCES_LPXX= $(TRAME_SRC_DIR)/lp/generic_lp.cpp
OBJECTS_LPXX= $(TRAME_SRC_DIR)/lp/generic_lp.o

# general flags
CFLAGS = -Wall $(OPT_FLAGS) $(LP_CFLAGS) -I$(TRAME_HEADER_DIR)
CXXFLAGS = $(CXX_STD) -Wall $(OPT_FLAGS) $(ARMA_INCLUDE_PATH) $(LP_CXXFLAGS) -I$(TRAME_HEADER_DIR)
FFLAGS = $(OPT_FLAGS)
LIBS= -lgfortran -framework Accelerate

# core TraME files
OBJECTS_TRAME_F90= $(TRAME_SRC_DIR)/prob/prob.o $(TRAME_SRC_DIR)/math/quadpack_double.o $(TRAME_SRC_DIR)/prob/trame_F90_aux.o

SOURCES_TRAME_AUX= $(TRAME_SRC_DIR)/aux/inv_pwa.cpp $(TRAME_SRC_DIR)/aux/logit_transform.cpp $(TRAME_SRC_DIR)/aux/trame_stats.cpp $(TRAME_SRC_DIR)/aux/zeroin.cpp
OBJECTS_TRAME_AUX= $(SOURCES_TRAME_AUX:.cpp=.o)

SOURCES_TRAME_ARUMS= $(TRAME_SRC_DIR)/arums/arums_empirical.cpp $(TRAME_SRC_DIR)/arums/arums_logit.cpp $(TRAME_SRC_DIR)/arums/arums_none.cpp $(TRAME_SRC_DIR)/arums/arums_probit.cpp $(TRAME_SRC_DIR)/arums/arums_rsc.cpp $(TRAME_SRC_DIR)/arums/arums_rusc.cpp
OBJECTS_TRAME_ARUMS= $(SOURCES_TRAME_ARUMS:.cpp=.o)

SOURCES_TRAME_MKTS= $(TRAME_SRC_DIR)/markets/dse.cpp $(TRAME_SRC_DIR)/markets/mmf.cpp $(TRAME_SRC_DIR)/markets/transfers.cpp
OBJECTS_TRAME_MKTS= $(SOURCES_TRAME_MKTS:.cpp=.o)

SOURCES_TRAME_MDLS= $(TRAME_SRC_DIR)/models/affinity.cpp $(TRAME_SRC_DIR)/models/model.cpp
OBJECTS_TRAME_MDLS= $(SOURCES_TRAME_MDLS:.cpp=.o)

SOURCES_TRAME_OPTIM= $(TRAME_SRC_DIR)/optim/bfgs.cpp $(TRAME_SRC_DIR)/optim/broyden.cpp $(TRAME_SRC_DIR)/optim/line_search.cpp $(TRAME_SRC_DIR)/optim/sumt.cpp $(TRAME_SRC_DIR)/optim/generic_optim.cpp $(TRAME_SRC_DIR)/optim/generic_constr_optim.cpp
OBJECTS_TRAME_OPTIM= $(SOURCES_TRAME_OPTIM:.cpp=.o)

SOURCES_TRAME_SLVRS= $(TRAME_SRC_DIR)/solvers/arc_newton.cpp $(TRAME_SRC_DIR)/solvers/aux_solvers.cpp $(TRAME_SRC_DIR)/solvers/eap_nash.cpp $(TRAME_SRC_DIR)/solvers/max_welfare.cpp $(TRAME_SRC_DIR)/solvers/nodal_newton.cpp
OBJECTS_TRAME_SLVRS= $(SOURCES_TRAME_SLVRS:.cpp=.o)

OBJECTS_TRAME= $(OBJECTS_LP) $(OBJECTS_LPXX) $(OBJECTS_TRAME_F90) $(OBJECTS_TRAME_AUX) $(OBJECTS_TRAME_ARUMS) $(OBJECTS_TRAME_MKTS) $(OBJECTS_TRAME_MDLS) $(OBJECTS_TRAME_OPTIM) $(OBJECTS_TRAME_SLVRS)

all: $(TRAME_DIR)/$(SHLIB) $(OBJECTS_TRAME)

# LP
$(OBJECTS_LP): $(SOURCES_LP)
	$(CC) $(CFLAGS) $^ -c -o $@

$(OBJECTS_LPXX): $(SOURCES_LPXX)
	$(CXX) $(CXXFLAGS) $^ -c -o $@

# aux and Fortran code
$(TRAME_SRC_DIR)/prob/prob.o: $(TRAME_SRC_DIR)/prob/prob.f90
	$(FC) $(FFLAGS) $^ -c -o $@

$(TRAME_SRC_DIR)/math/quadpack_double.o: $(TRAME_SRC_DIR)/math/quadpack_double.f90
	$(FC) $(FFLAGS) $^ -c -o $@

$(TRAME_SRC_DIR)/prob/trame_F90_aux.o: $(TRAME_SRC_DIR)/prob/trame_F90_aux.f90
	$(FC) $(FFLAGS) $^ -c -o $@

$(TRAME_SRC_DIR)/aux/%.o: $(TRAME_SRC_DIR)/aux/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(TRAME_SRC_DIR)/optim/%.o: $(TRAME_SRC_DIR)/optim/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

# core TraME files
$(TRAME_SRC_DIR)/arums/%.o: $(TRAME_SRC_DIR)/arums/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(TRAME_SRC_DIR)/markets/%.o: $(TRAME_SRC_DIR)/markets/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(TRAME_SRC_DIR)/models/%.o: $(TRAME_SRC_DIR)/models/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

$(TRAME_SRC_DIR)/solvers/%.o: $(TRAME_SRC_DIR)/solvers/%.cpp
	$(CXX) $(CXXFLAGS) $< -c -o $@

# shared library
$(TRAME_DIR)/$(SHLIB): $(OBJECTS_TRAME)
	$(CXX) $(SHLIB_FLAGS) -o $@ $^ $(LP_LIBS) $(LIBS)
#    ar -rvs ../../inst/libs/libtrame.a $^

# cleanup and install
.PHONY: clean
clean:
	@rm -f *.so ./tests/*/*.test ./tests/*/*.o $(TRAME_SRC_DIR)/*/*.o

.PHONY: install
install: $(SHLIB)
	@cp $(TRAME_DIR)/$(SHLIB) $(INSTALL_PATH)/lib/$(SHLIB)
	@mkdir -p $(INSTALL_PATH)/include/trame
	@cp -r $(TRAME_DIR)/include/* $(INSTALL_PATH)/include/trame
