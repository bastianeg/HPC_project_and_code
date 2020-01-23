#ifndef __TRANSFER_3D_GPU
#define __TRANSFER_3D_GPU

void transfer_3d(double ***dst, double ***src, int nx, int ny, int nz, int flag);
void transfer_3d_to_1d(double *dst, double ***src, int nx, int ny, int nz, int flag);
void transfer_3d_from_1d(double ***dst, double *src, int nx, int ny, int nz, int flag);

#endif /* __TRANSFER_3D_GPU */
