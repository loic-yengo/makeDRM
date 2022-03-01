makeDRM: 	makeDRM.cpp
	g++ -fopenmp -I"./eigen-3.4.0/" -Wall -std=c++11 -g -O2 makeDRM.cpp -o makeDRM -lm

makeDRMfloat:        makeDRM_float.cpp
	g++ -fopenmp -I"./eigen-3.4.0/" -Wall -std=c++11 -g -O2 makeDRM_float.cpp -o makeDRMfloat -lm

clean:
	rm makeDRM; rm makeDRMfloat

