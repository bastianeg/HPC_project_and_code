TARGET	= libmatmult.so
LIBSRCS	= matmult_blk.c native_combined.c
LIBOBJS	= matmult_blk.o native_combined.o

OPT	= -g -O3 -ffast-math -funroll-loops
PIC	= -fPIC

CC	= gcc
CFLAGS= $(OPT) $(PIC) $(XOPTS)

SOFLAGS = -shared 
XLIBS	= 

$(TARGET): $(LIBOBJS)
	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS)

clean:
	@/bin/rm -f core core.* $(LIBOBJS) 
