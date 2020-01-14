CC = gcc
CFLAGS = -Wall 
LDLIBS  = -lopenblas
LDFLAGS = -L /usr/lib64/atlas -lsatlas
OBJECTS = matrixmatrix.o
excer1: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o matrixmatrix $(OBJECTS)

clean:
	@rm -f *.o core