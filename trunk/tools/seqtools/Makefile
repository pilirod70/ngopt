all: fastalen revcomp rmctgs rmfanl exmap

fastalen: fastalen.cpp
	g++ -o fastalen fastalen.cpp

revcomp: revcomp.cpp
	g++ -o revcomp revcomp.cpp

rmctgs: rmctgs.cpp
	g++ -o rmctgs rmctgs.cpp

rmfanl: rmfastanl.cpp
	g++ -o rmfanl rmfastanl.cpp

exmap: extract-mapped.cpp
	g++ -O3 -static -std=c++0x -o exmap  extract-mapped.cpp
