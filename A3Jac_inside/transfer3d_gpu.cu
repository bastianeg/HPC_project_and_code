#include <cuda_runtime_api.h>
#include <helper_cuda.h>

extern "C" {

    void
    transfer_3d(double ***dst, double ***src, int nx, int ny, int nz, int flag)
    {
        long nPtr = nz + nz * ny;
        long nBlk = nx * ny * nz;

        // we only transfer the value block
        checkCudaErrors( cudaMemcpy((double *) dst + nPtr,
                                    (double *) src + nPtr,
                                    nBlk * sizeof(double),
                                    (cudaMemcpyKind) flag) );
    }

    void
    transfer_3d_to_1d(double *dst, double ***src, int nx, int ny, int nz, int flag)
    {
        long nPtr = nz + nz * ny;
        long nBlk = nx * ny * nz;

        // we only transfer the value block
        checkCudaErrors( cudaMemcpy((double *) dst,
                                    (double *) src + nPtr,
                                    nBlk * sizeof(double),
                                    (cudaMemcpyKind) flag) );
    }

    void
    transfer_3d_from_1d(double ***dst, double *src, int nx, int ny, int nz, int flag)
    {
        long nPtr = nz + nz * ny;
        long nBlk = nx * ny * nz;

        // we only transfer the value block
        checkCudaErrors( cudaMemcpy((double *) dst + nPtr,
                                    (double *) src,
                                    nBlk * sizeof(double),
                                    (cudaMemcpyKind) flag) );
    }

}
