# Makefile
# hello
TARGET_S  = poisson_seq		# Jacobi
TARGET_N = poisson_naive		# Gauss-Seidel
TARGET_MUL = poisson_multi		# Jacobi parallel
TARGET_TOL = poisson_tol		# Gauss-Seidel parallel

SOURCES	= main.c print.c alloc3d.c
OBJECTS	= print.o alloc3d.o
MAIN_S	= main_seq.o
MAIN_N = main_n.o
MAIN_MUL= main_mul.o
MAIN_TOL = main_tol.o
OBJS_S	= $(MAIN_S) jacobiseq.o
OBJS_N	= $(MAIN_N) jacobinaive.o
OBJS_M	= $(MAIN_MUL) jacobimulti.o
OBJS_TOL	= $(MAIN_TOL) jacobitol.o

# options and settings for the GCC compilers
#
CXX	= nvcc
DEFS	=
OPT	= -g -O3 -funroll-loops -ffast-math
IPO	=
ISA	=
CHIP	=
ARCH	=
PARA	=
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS) -fopenmp
LDFLAGS = -lm

all: $(TARGET_S) $(TARGET_N) $(TARGET_MUL) $(TARGET_TOL)

$(TARGET_S): $(OBJECTS) $(OBJS_S)
	$(CC) -o $@ $(CFLAGS) $(OBJS_S) $(OBJECTS) $(LDFLAGS)

$(TARGET_N): $(OBJECTS) $(OBJS_N)
	$(CC) -o $@ $(CFLAGS) $(OBJS_N) $(OBJECTS) $(LDFLAGS)

$(TARGET_MUL): $(OBJECTS) $(OBJS_MUL)
	$(CC) -o $@ $(CFLAGS) $(OBJS_MUL) $(OBJECTS) $(LDFLAGS)

$(TARGET_TOL): $(OBJECTS) $(OBJS_TOL)
	clang -o $@ $(CFLAGS) $(OBJS_TOL) $(OBJECTS) $(LDFLAGS)

$(MAIN_S):
	$(CC) -o $@ -D_JACOBI $(CFLAGS) -c main.c

$(MAIN_N):
	$(CC) -o $@ -D_GAUSS_SEIDEL $(CFLAGS) -c main.c

$(MAIN_MUL):
	$(CC) -o $@ -D_JACOBI_PAR $(CFLAGS) -c main.c

$(MAIN_TOL):
	$(CC) -o $@ -D_GAUSS_SEIDEL_PAR $(CFLAGS) -c main.c

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET_S)  $(TARGET_N) $(TARGET_MUL)  $(TARGET_TOL)

# DO NOT DELETE

main_j.o: main.c print.h jacobiseq.h
main_n.o: main.c print.h jacobinaive.h
main_mul.o: main.c print.h jacobimulti.h
main_tol.o: main.c print.h jacobitol.h
print.o: print.h
