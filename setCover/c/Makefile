CFLAGS= -O0 -Wall -Wextra -std=c99 -mtune=native -fopenmp
CC=gcc

MODE?=NONE

set_cover: set_cover.o protein_oligo_library.o dynamic_string.o hash_table.o array_list.o set.o 
	$(CC) $(CFLAGS) -D$(MODE) -o $@ $^
set_cover.o: set_cover.c protein_oligo_library.h hash_table.h array_list.h set.h
	$(CC) $(CFLAGS) -D$(MODE) -c -o $@ $<

protein_oligo_library.o: protein_oligo_library.c protein_oligo_library.h hash_table.h array_list.h set.h
	$(CC) $(CFLAGS) -c -o $@ $<

dynamic_string.o: dynamic_string.c dynamic_string.h
	$(CC) $(CFLAGS) -c -o $@ $<

hash_table.o: hash_table.c hash_table.h
	$(CC) $(CFLAGS) -c -o $@ $<

array_list: array_list_main.o array_list.o
	$(CC) $(CFLAGS) -c -o $@ $<
array_list_main.o: array_list_main.c array_list.h
	$(CC) $(CFLAGS) -c -o $@ $<
array_list.o: array_list.c array_list.h
	$(CC) $(CFLAGS) -c -o $@ $<

set.o: set.c set.h 
	$(CC) $(CFLAGS) -c -o $@ $<


.PHONY: debug clean optimized profile time_trial
debug: CFLAGS+= -g -O0
debug: clean
debug: set_cover

optimized: CFLAGS += -O3  -ffast-math
optimized: clean
optimized: set_cover

time_trial: clean
time_trial: MODE=TIME_TRIAL
time_trial: optimized

profile: CFLAGS += -pg -g
profile: clean
profile: set_cover


clean:
	rm -rf *.o *.gch set_cover
