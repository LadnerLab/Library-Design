CFLAGS= -O0 -Wall -fopenmp 
CC= g++

main: main.o 
	$(CC) $(CFLAGS) main.o -o main
main.o: main.cpp xoroshiro.h table.h
	$(CC) $(CFLAGS) main.cpp -c -o main.o


.PHONY: clean optimized debug profile

clean: 
	rm -f *.o main *.gch

optimized: CFLAGS+= -O3 -mtune=native 
optimized: clean
optimized: main

profile: CFLAGS+= -pg 
profile: clean
profile: main

debug: clean
debug: CFLAGS+= -g
debug: main
