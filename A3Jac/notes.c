_global_ void
setreszero(int n, double* result){
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    if(tid<n){
        result[tid] = 0.0;
    }
}


_global_ void 
mvgpu(int n, double* A, double* b, double* result){

    int i = threadIdx.x

    for(int i = 0 ; i < n; i++){
        for(int j = 0; j < n; j++){
            result[i] += b[j]*A[j+n*i];
        }
    }
}

setreszero<<<N/BLOCK_DIM,BLOCK_DIM>>>(N, d_result);
        cudaDeviceSynchronize();
        mvgpu<<<dim3 (N/BLOCK_DIM,N/BLOCK_DIM),dim3 (BLOCK_DIM,BLOCK_DIM)>>>(N,d_A,d_b,d_result);
        cudaDeviceSynchronize();

// allocate memeory for image on host
    cudaMallocHost((void **) &h_image, IM_SIZE * IM_SIZE * sizeof(int));
    
    if ( h_image == NULL ) {
       fprintf(stderr, "memory allocation failed!\n");
       return(1);
    }
    // allocate memory for image on device
    cudaMalloc((void **) &d_image, IM_SIZE * IM_SIZE * sizeof(int));
    // move to device
    double t0 = omp_get_wtime();
    cudaMemcpy(d_image,h_image,IM_SIZE * IM_SIZE * sizeof(int),cudaMemcpyHostToDevice);

    // start kernel
    mandelgpu<<<dim3 (IM_SIZE/BLOCK_DIM,IM_SIZE/BLOCK_DIM),dim3 (BLOCK_DIM,BLOCK_DIM)>>>(IM_SIZE, IM_SIZE, d_image, max_iter);
    cudaDeviceSynchronize();
    
    // copy back to host
    cudaMemcpy(h_image,d_image,IM_SIZE * IM_SIZE * sizeof(int),cudaMemcpyDeviceToHost);
    double elapsed = omp_get_wtime() - t0;
    // print image
    printf("time elapsed %.8f\n",elapsed);
    writepng("mandelbrot.png", h_image, IM_SIZE, IM_SIZE);

    // free memory on both host and device
    cudaFreeHost(h_image); cudaFree(d_image);

    _global_ void
mandelgpu(int disp_width, int disp_height, int *array, int max_iter) {

	int j = blockIdx.x*blockDim.x+threadIdx.x;
	int i = blockIdx.y*blockDim.y+threadIdx.y;

	if((j<disp_width)&&(i<disp_height)){
		double scale_real = 3.5 / (double)disp_width;
		double scale_imag = 3.5 / (double)disp_height;

		double x = ((double)i * scale_real) - 2.25; 
		double y = ((double)j * scale_imag) - 1.75; 

		double u    = 0.0;
		double v    = 0.0;
		double u2   = 0.0;
		double v2   = 0.0;
		int iter = 0;

		while ( u2 + v2 < 4.0 &&  iter < max_iter ) {
			v = 2 * v * u + y;
			u = u2 - v2 + x;
			u2 = u*u;
			v2 = v*v;
			iter = iter + 1;
		}

		// if we exceed max_iter, reset to zero
		iter = iter == max_iter ? 0 : iter;

		array[i*disp_height + j] = iter;
	}
}