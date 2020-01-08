TARGET	= libmatmult.so
LIBSRCS	= matmult_blk.c native_combined.c cblas.c
LIBOBJS	= matmult_blk.o native_combined.o cblas.o
LDFLAGS = -L /usr/lib64/atlas -lsatlas

OPT	= -g -O3 -ffast-math -funroll-loops
PIC	= -fPIC

CC	= gcc
CFLAGS= -std=c99 $(OPT) $(PIC) $(XOPTS) 

SOFLAGS = -shared 
XLIBS	= 

$(TARGET): $(LIBOBJS)
	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS)

clean:
	@/bin/rm -f core core.* $(LIBOBJS) 
