CXX = g++
LD = $(CXX)

OPTIM = -O0
CXXFLAGS += -Wall -pedantic -std=c++17 $(OPTIM)

INCLUDE = -I ../../include
LIBPATH = ../../lib
LIBNAME = SHTOOLS
FFTW = -lfftw3 -lm
LAPACK = -llapack
BLAS = -lblas

LIB = -L $(LIBPATH) -l$(LIBNAME) -lgfortran -ldl $(FFTW) $(LAPACK) $(BLAS)

all : test

test.o : test.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<
	
test : test.o
	$(LD) -o $@ $< $(LIB)
	
run-cpp-tests-no-timing :
	./test
	
clean :
	@-rm -rf *.o
	@-rm -rf test
