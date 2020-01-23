# Makefile
# hello
TARGET_S  = jacobiseq		# Jacobi
TARGET_N = jacobinaive		# Naive
TARGET_MUL = jacobimulti	# Multi
TARGET_TOL = jacobitol		# Tolerance

SOURCES	= main.cu print.c
OBJECTS	= print.o
MAIN_S	= main_seq.o
MAIN_N = main_n.o
MAIN_MUL= main_mul.o
MAIN_TOL = main_tol.o
OBJS_S	= $(MAIN_S) jacobiseq.o
OBJS_N	= $(MAIN_N) jacobinaive.o
OBJS_M	= $(MAIN_MUL) jacobimulti.o
OBJS_TOL= $(MAIN_TOL) jacobitol.o

# options and settings for the GCC compilers
#

PIC = -fpic
OMP   = -fopenmp
XPIC  = -Xcompiler -fpic
XOPT  = -Xptxas=-v # use -lineinfo for profiler, use -G for debugging
XARCH = -arch=sm_70


CXX	= nvcc
OPT	= -g -O3
CXXFLAGS	= --compiler-options "-fopenmp" $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS)

CUDA_PATH ?= /appl/cuda/10.2
INCLUDES = -I$(CUDA_PATH)/include \
	       -I$(CUDA_PATH)/samples/NVIDIA_CUDA-10.2_Samples/common/inc

LDFLAGS = -lm

.SUFFIXES: .cu
.cu.o:
	$(CXX) -o $*.o -c $*.cu $(CXXFLAGS) $(SOFLAGS) $(INCLUDES)

all: $(TARGET_S) $(TARGET_N) $(TARGET_MUL) $(TARGET_TOL)

$(TARGET_S): $(OBJECTS) $(OBJS_S)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS_S) $(OBJECTS) $(LDFLAGS)

$(TARGET_N): $(OBJECTS) $(OBJS_N)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS_N) $(OBJECTS) $(LDFLAGS)

$(TARGET_MUL): $(OBJECTS) $(OBJS_MUL)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS_MUL) $(OBJECTS) $(LDFLAGS)

$(TARGET_TOL): $(OBJECTS) $(OBJS_TOL)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS_TOL) $(OBJECTS) $(LDFLAGS)
 
$(MAIN_S):
	$(CXX) -o $@ -D_JACOBISEQ $(CXXFLAGS) -c main.cu

$(MAIN_N):
	$(CXX) -o $@ -D_JACOBINAIVE $(CXXFLAGS) -c main.cu

$(MAIN_MUL):
	$(CXX) -o $@ -D_JACOBIMULTI $(CXXFLAGS) -c main.cu

$(MAIN_TOL):
	$(CXX) -o $@ -D_JACOBITOL $(CXXFLAGS) -c main.cu

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET_S)  $(TARGET_N) $(TARGET_MUL)  $(TARGET_TOL)

# DO NOT DELETE

main_j.o: main.cu print.h jacobiseq.h
main_n.o: main.cu print.h jacobinaive.h
main_mul.o: main.cu print.h jacobimulti.h
main_tol.o: main.cu print.h jacobitol.h
print.o: print.h
