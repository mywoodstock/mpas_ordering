
NETCDF = /usr
#BOOST = /usr/local/boost

CXX = g++

OPT_FLAGS = -O3
DEBUG_FLAGS = -g -DDEBUG
CXX_FLAGS = -std=c++11 -D_64BITOFFSET

CXX_INCL = -I$(NETCDF)/include #-I$(BOOST)/include
CXX_LIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf

OBJ = netcdf_utils.o mpas_ordering.o mpas_order.o mpas_ordering_utils.o \
	  mpas_ordering_xyzsort.o mpas_ordering_random.o mpas_ordering_morton.o \
	  mpas_ordering_hilbert.o mpas_ordering_peano.o

BIN = mpas_order

DEPS_HPP = mpas_ordering.hpp mpas_ordering_typedefs.hpp netcdf_utils.h

all: CXX_FLAGS+=$(OPT_FLAGS)
all: $(BIN)

debug: CXX_FLAGS+=$(DEBUG_FLAGS)
debug: $(BIN)

clean:
	rm -f $(OBJ) $(BIN)


$(BIN): $(OBJ)
	$(CXX) -o $@ $^ $(CXX_FLAGS) $(CXX_LIBS)

%.o: %.cpp $(DEPS_HPP)
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)


#netcdf_utils.o: netcdf_utils.cpp netcdf_utils.h
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
#
#mpas_ordering.o: mpas_ordering.cpp mpas_ordering.hpp
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
#
#mpas_ordering_xyzsort.o: mpas_ordering_xyzsort.cpp mpas_ordering.hpp
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
#
#mpas_ordering_random.o: mpas_ordering_random.cpp mpas_ordering.hpp
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
#
#mpas_ordering_morton.o: mpas_ordering_morton.cpp mpas_ordering.hpp
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
#
#mpas_order.o: mpas_order.cpp
#	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
