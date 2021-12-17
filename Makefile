.PHONY: all

MPICXX := mpicxx
CXXFLAGS += -g -std=c++11 -O3 -ffast-math -march=native

OBJS := BOBHash32.o mrac.o mpi-mrac-v1.o mpi-mrac-v2.o
EXES := mrac mpi-mrac-v1 mpi-mrac-v2

all: $(EXES)

mrac: BOBHash32.o mrac.o
	$(CXX) $(CXXFLAGS) $^ -o $@

mpi-mrac-v1: BOBHash32.o mpi-mrac-v1.o
	$(MPICXX) $(CXXFLAGS) $^ -o $@

mpi-mrac-v2: BOBHash32.o mpi-mrac-v2.o
	$(MPICXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f $(OBJS)
	rm -f $(EXES)