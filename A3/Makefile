TARGET	= libmatmult.so
LIBSRCS	= matmult_gpun.cu matmult_kernels.cu assign1_functions.cu
LIBOBJS	= matmult_gpun.o matmult_kernels.o assign1_functions.o

OPT	= -g -O3
PIC = -fpic
OMP   = -fopenmp
XPIC  = -Xcompiler -fpic
XOPT  = -Xptxas=-v # use -lineinfo for profiler, use -G for debugging
XARCH = -arch=sm_70

CXX	= nvcc
CXXFLAGS = --compiler-options "$(OPT) $(PIC) $(OMP)" $(XARCH) $(XOPT) $(XPIC) 

CUDA_PATH ?= /appl/cuda/10.2
INCLUDES = -I$(CUDA_PATH)/include \
	       -I$(CUDA_PATH)/samples/NVIDIA_CUDA-10.2_Samples/common/inc

SOFLAGS = -shared
XLIBS	= -lcublas

$(TARGET): $(LIBOBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(SOFLAGS) $(INCLUDES) $^ $(XLIBS)

.SUFFIXES: .cu
.cu.o:
	$(CXX) -o $*.o -c $*.cu $(CXXFLAGS) $(SOFLAGS) $(INCLUDES)

clean:
	/bin/rm -f $(TARGET) $(LIBOBJS) 

matmult_gpun.o: matmult_gpun.cu matmult_kernels.h
matmult_kernels.o: matmult_kernels.cu
assign1_functions.o: assign1_functions.cu
