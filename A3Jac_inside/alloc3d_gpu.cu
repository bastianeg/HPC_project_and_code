#include <cuda_runtime_api.h>
#include <helper_cuda.h>

__global__ void d_malloc_3d_gpu_kernel1(double *** array3D, int nx, int ny, int nz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nz) {
        array3D[i] = (double **) array3D + nz + i * ny;
        //printf("k1: %i | %i\n", i, i* ny);
    }
}

__global__ void d_malloc_3d_gpu_kernel2(double *** array3D, int nx, int ny, int nz)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < nz && j < ny) {
        array3D[i][j] = (double *) array3D + nz + nz * ny + (i * nx * ny) + (j * nx);
        //printf("k2: %i %i | %i\n", i, j, (i * nx * ny) + (j * nx));
    }
}

extern "C" {

    double ***
    d_malloc_3d_gpu(int nx, int ny, int nz) {

        if (nx <= 0 || ny <= 0 || nz <= 0)
            return NULL;
    
        double ***array3D; 
        checkCudaErrors( cudaMalloc((void**)&array3D, 
                                    nz * sizeof(double **) +
                                    nz * ny * sizeof(double *) +
                                    nz * ny * nx * sizeof(double)) );
        if (array3D == NULL) {
            return NULL;
        }

        dim3 block(16, 16);
        dim3 grid((nz + 15) / 16, (ny + 15) /16);
        d_malloc_3d_gpu_kernel1<<<grid.x, block.x>>>(array3D, nx, ny, nz);
        d_malloc_3d_gpu_kernel2<<<grid, block>>>(array3D, nx, ny, nz);
        checkCudaErrors( cudaDeviceSynchronize() );

        return array3D;
    }

    void
    free_gpu(double ***array3D) {
        checkCudaErrors( cudaFree(array3D) );
    }

}
