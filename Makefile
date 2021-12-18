.PHONY: all

MPICXX := mpicxx
CXXFLAGS += -g -std=c++11 -O3 -ffast-math -march=native

OBJS := BOBHash32.o mrac.o mrac-buf.o mpi-mrac-v1.o mpi-mrac-v2.o
EXES := mrac mrac-buf mpi-mrac-v1 mpi-mrac-v2

all: $(EXES)

mrac: BOBHash32.o mrac.o
	$(CXX) $(CXXFLAGS) $^ -o $@

mrac-buf: BOBHash32.o mrac-buf.o
	$(CXX) $(CXXFLAGS) $^ -o $@

mpi-mrac-v1: BOBHash32.o mpi-mrac-v1.o
	$(MPICXX) $(CXXFLAGS) $^ -o $@

mpi-mrac-v2: BOBHash32.o mpi-mrac-v2.o
	$(MPICXX) $(CXXFLAGS) $^ -o $@

mrac.o: mrac.cpp EMFSD.h MRAC.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

mrac-buf.o: mrac.cpp EMFSD.h MRAC.h
	$(CXX) $(CXXFLAGS) -DUSE_BUFFER $< -c -o $@

mpi-mrac-v1.o: mpi-mrac-v1.cpp EMFSD.h MRAC.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

mpi-mrac-v2.o: mpi-mrac-v2.cpp EMFSD.h MRAC.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

BOBHash32.o: BOBHash32.cpp BOBHash32.h
	$(CXX) $(CXXFLAGS) $< -c -o $@

clean:
	rm -f $(OBJS)
	rm -f $(EXES)