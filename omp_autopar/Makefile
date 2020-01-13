CC	= gcc
OPT	= -g -O3 -ffast-math
CHIP	= 
ISA	= 
PARALELL= -ftree-parallelize-loops=4
LIBS	= -lm -xopenmp
CFLAGS  = $(OPT) $(CHIP) $(ISA) $(PARALELL) $(XOPT)

SOURCES = calc_pi.c 
OBJECTS = $(SOURCES:.c=.o)

TARGET = progam.$(CC)

all:
	$(MAKE) $(MAKEFLAGS) clean
	$(MAKE) $(MAKEFLAGS) $(TARGET)

$(TARGET): clean $(OBJECTS)
	$(CC) -o $@ $(CFLAGS) $(OBJECTS) $(LIBS)

clean:
	@/bin/rm -f $(OBJECTS) core

realclean: clean
	@/bin/rm -f $(TARGET)
