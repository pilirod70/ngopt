all: gc_calc merge_gc-plp calc_avg-gc calc_avg-cov

gc_calc: gc_calc.cpp
	g++ -o gc_calc gc_calc.cpp

merg_gc-plp: merge_gc-plp.cpp
	g++ -o merge_gc-plp merge_gc-plp.cpp

calc_avg-gc: calc_avg-gc.cpp
	g++ -o calc_avg-gc calc_avg-gc.cpp

calc_avg-cov: calc_avg-cov.cpp
	g++ -o calc_avg-cov calc_avg-cov.cpp

clean:
	rm gc_calc merge_gc-plp calc_avg-gc calc_avg-cov
