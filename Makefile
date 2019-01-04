CPP=g++ -std=c++11
CPPFLAGS=-O3 -Wall -DVERBOSE
SOURCES=./extpoints.cpp

unary:
	@echo "Extreme Points Search for Convex Hull with openMP improvements"
	@$(CPP) -fopenmp $(CPPFLAGS) $(INCLUDES) $(SOURCES) -o prog
