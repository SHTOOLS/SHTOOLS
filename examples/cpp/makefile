CXX = g++
LD = $(CXX)

OPTIM = -O0
CXXFLAGS += -Wall -pedantic -std=c++17 $(OPTIM) 

INCLUDE=-I ../../src

LIB= -L ../../lib -lSHTOOLS -lgfortran \
     -ldl -lfftw3 -lm -llapack -lblas \

all: test

test.o : test.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $<
	
test : test.o 
	$(LD) -o $@ $< $(LIB)
	
clean :
	rm -rf *.o
	rm -rf test
