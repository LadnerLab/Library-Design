CFLAGS= -O0 -Wall -Wextra -std=c99 -fopenmp

one_hundred_reps: one_hundred_reps.o array_list.o dynamic_string.o protein_oligo_library.o
	gcc $(CFLAGS) one_hundred_reps.o array_list.o dynamic_string.o protein_oligo_library.o -lpthread -o one_hundred_reps
one_hundred_reps.o: one_hundred_reps.c array_list.h dynamic_string.h protein_oligo_library.h
array_list.o: array_list.c array_list.h
dynamic_string.o: dynamic_string.c dynamic_string.h
protein_oligo_library.o: protein_oligo_library.c protein_oligo_library.h

.PHONY: debug optimized
debug: CFLAGS += -g 
debug: clean
debug: one_hundred_reps

optimized: CFLAGS += -O3
optimized: clean
optimized: one_hundred_reps


clean:
	rm -f *.o *.gch one_hundred_reps *.~
