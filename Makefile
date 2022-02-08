makeDRM: 	makeDRM.cpp
	g++ -fopenmp -I"./eigen-3.4.0/" -Wall -std=c++11 -g -O2 makeDRM.cpp -o makeDRM -lm

clean:
	rm makeDRM

