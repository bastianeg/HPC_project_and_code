CC	= suncc
OPT	= -g -fast -xrestrict
CHIP	=
ISA	=
PARALELL= #-xautopar -xloopinfo -xreduction
LIBS	= -lm -xopenmp
CFLAGS  = $(OPT) $(CHIP) $(ISA) $(PARALELL) $(XOPT) -xopenmp

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
