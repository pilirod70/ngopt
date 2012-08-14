all: pairfa pairfq repair

pairfa: pairfa.cpp 
	g++ -O3 -static -std=c++0x -o pairfa pairfa.cpp

pairfq: pairfq.cpp
	g++ -O3 -static -std=c++0x -o pairfq pairfq.cpp

repair: repair.cpp
	g++ -O3 -static -std=c++0x -o repair repair.cpp

clean: 
	rm pairfa pairfq repair
#g++ -O3 -static -std=c++0x -o pairfq pairfq.cpp

