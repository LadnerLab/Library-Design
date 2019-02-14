CFLAGS= -g -O0 
CC= g++

main: main.o 
	$(CC) $(CFLAGS) main.o -o main
main.o: main.cpp xoroshiro.h table.h
	$(CC) $(CFLAGS) main.cpp -c -o main.o


clean:
	rm -f *.o main
