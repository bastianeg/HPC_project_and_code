CC = gcc
OPT = -g -O3 -ffast-math -funroll-loops
CFLAGS = -Wall 
LDFLAGS = 
OBJECTS = test.o matmult_blk.o
excer1: $(OBJECTS)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o test $(OBJECTS)

clean:
	@rm -f *.o core

matmul_blk.o: matmult_blk.h