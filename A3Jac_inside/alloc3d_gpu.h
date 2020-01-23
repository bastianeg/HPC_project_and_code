#ifndef __ALLOC_3D_GPU
#define __ALLOC_3D_GPU

double ***d_malloc_3d_gpu(int m, int n, int k);
void free_gpu(double ***array3D);

#endif /* __ALLOC_3D_GPU */
