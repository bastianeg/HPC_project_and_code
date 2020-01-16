# Makefile
# hello
TARGET_J  = poisson_j		# Jacobi
TARGET_GS = poisson_gs		# Gauss-Seidel
TARGET_JP  = poisson_jp		# Jacobi parallel
TARGET_GSP = poisson_gsp		# Gauss-Seidel parallel

SOURCES	= main.c print.c alloc3d.c
OBJECTS	= print.o alloc3d.o
MAIN_J	= main_j.o
MAIN_GS = main_gs.o
MAIN_JP	= main_jp.o
MAIN_GSP = main_gsp.o
OBJS_J	= $(MAIN_J) jacobi.o
OBJS_GS	= $(MAIN_GS) gauss_seidel.o
OBJS_JP	= $(MAIN_JP) jacobi_par.o
OBJS_GSP	= $(MAIN_GSP) gauss_seidel_par.o

# options and settings for the GCC compilers
#
CC	= clang
DEFS	=
OPT	= -g -O3 -funroll-loops -ffast-math
IPO	=
ISA	=
CHIP	=
ARCH	=
PARA	=
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS) -fopenmp
LDFLAGS = -lm

all: $(TARGET_J) $(TARGET_GS) $(TARGET_JP) $(TARGET_GSP)

$(TARGET_J): $(OBJECTS) $(OBJS_J)
	$(CC) -o $@ $(CFLAGS) $(OBJS_J) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS): $(OBJECTS) $(OBJS_GS)
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS) $(OBJECTS) $(LDFLAGS)

$(TARGET_JP): $(OBJECTS) $(OBJS_JP)
	$(CC) -o $@ $(CFLAGS) $(OBJS_JP) $(OBJECTS) $(LDFLAGS)

$(TARGET_GSP): $(OBJECTS) $(OBJS_GSP)
	$(CC) -o $@ $(CFLAGS) $(OBJS_GSP) $(OBJECTS) $(LDFLAGS)

$(MAIN_J):
	$(CC) -o $@ -D_JACOBI $(CFLAGS) -c main.c

$(MAIN_GS):
	$(CC) -o $@ -D_GAUSS_SEIDEL $(CFLAGS) -c main.c

$(MAIN_JP):
	$(CC) -o $@ -D_JACOBI_PAR $(CFLAGS) -c main.c

$(MAIN_GSP):
	$(CC) -o $@ -D_GAUSS_SEIDEL_PAR $(CFLAGS) -c main.c

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET_J)  $(TARGET_GS) $(TARGET_JP)  $(TARGET_GSP)

# DO NOT DELETE

main_j.o: main.c print.h jacobi.h
main_gs.o: main.c print.h gauss_seidel.h
main_jp.o: main.c print.h jacobi_par.h
main_gsp.o: main.c print.h gauss_seidel_par.h
print.o: print.h
