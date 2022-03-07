CFLAGS= -O0 -Wall -fopenmp
CC= g++

oligo_encoding: main.o
	$(CC) $(CFLAGS) main.o -o main
main.o: main.cpp xoroshiro.h table.h
	$(CC) $(CFLAGS) main.cpp -c -o main.o


.PHONY: clean optimized debug profile

clean:
	rm -f *.o oligo_ecoding *.gch

optimized: CFLAGS+= -O3 -mtune=native
optimized: clean
optimized: oligo_ecoding

profile: CFLAGS+= -pg
profile: clean
profile: oligo_ecoding

debug: clean
debug: CFLAGS+= -g
debug: oligo_ecoding
