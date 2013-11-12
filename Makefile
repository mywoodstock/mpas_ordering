
NETCDF = /usr

CXX = g++

OPT_FLAGS = -O3
DEBUG_FLAGS = -g

CXX_INCL = -I$(NETCDF)/include
CXX_LIBS = -L$(NETCDF)/lib -lnetcdf_c++ -lnetcdf

OBJ = netcdf_utils.o mpas_ordering.o mpas_order.o

all:
	CXX_FLAGS = $(OPT_FLAGS)
	mpas_order

debug:
	CXX_FLAGS = $(DEBUG_FLAGS)
	mpas_order

clean:
	rm *.o mpas_order

mpas_order: $(OBJ)
	$(CXX) -o $@ $^ $(CXX_FLAGS) $(CXX_LIBS)

netcdf_utils.o: netcdf_utils.cpp netcdf_utils.h
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)

mpas_ordering.o: mpas_ordering.cpp mpas_ordering.hpp
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)

mpas_order.o: mpas_order.cpp
	$(CXX) -c $< -o $@ $(CXX_FLAGS) $(CXX_INCL)
