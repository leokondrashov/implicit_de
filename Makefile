CXX = g++ 
CXXFLAGS = -Wall -g -MMD

.PHONY: clean all

all: implicit 

clean:
	rm *.d
	rm *.o
	rm implicit

implicit: vector.o matrix.o

-include *.d

