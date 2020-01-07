CC = gcc
CFLAGS = -Wall 
LDFLAGS = -framework Accelerate
OBJECTS = matrixmatrix.o
excer1: $(OBJECTS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o excer1 $(OBJECTS)

clean:
	@rm -f *.o core