CC = g++ -ggdb



all : realall
	make bin

realall :
	(cd ./src; $(CC) -c point.cpp; $(CC) fish.cpp point.o -o fish)

bin : 
	mkdir bin
	mv src/fish ./bin
clean :
	(cd src ; rm *.o)
	(rm -rf ./bin)
