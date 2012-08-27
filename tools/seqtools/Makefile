all: spscaf fastalen revcomp rmctgs rmfanl rmfqnl exmap sumfasta

spscaf: split_scaffolds.cpp
	g++ -o spscaf split_scaffolds.cpp

fastalen: fastalen.cpp
	g++ -o fastalen fastalen.cpp

revcomp: revcomp.cpp
	g++ -o revcomp revcomp.cpp

rmctgs: rmctgs.cpp
	g++ -o rmctgs rmctgs.cpp

rmfanl: rmfastanl.cpp
	g++ -o rmfanl rmfastanl.cpp

rmfqnl: rmfastqnl.cpp
	g++ -o rmfqnl rmfastqnl.cpp

exmap: extract-mapped.cpp
	g++ -o exmap  extract-mapped.cpp

sumfasta: sumfasta.cpp
	g++ -o sumfasta sumfasta.cpp
