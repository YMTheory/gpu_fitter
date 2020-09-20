#include "cuda_runtime.h"
#include <stdlib.h>
#include <iostream>
#include <sys/time.h>

using namespace std;

__global__ void add(float* x, float* y, float* z, int n)
{
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int stride = blockDim.x * gridDim.x;
    for(int i=index; i<n; i+=stride) {
        z[i] = x[i] + y[i];
    }
}

int main()
{
    struct timeval start, end;
    gettimeofday( &start, NULL );

    int N = 1000; //1 << 20;
    int nBytes = N * sizeof(float);
    
    // CPU memory allocation
    float *x, *y, *z;
    x = (float*)malloc(nBytes);
    y = (float*)malloc(nBytes);
    z = (float*)malloc(nBytes);

    // array initialization
    for(int i=0; i<N; ++i) {
        x[i] = 10.0;
        y[i] = 20.0;
    }

    // GPU memory initialization
    float *dx, *dy, *dz ;
    cudaMalloc((void**)&dx, nBytes);
    cudaMalloc((void**)&dy, nBytes);
    cudaMalloc((void**)&dz, nBytes);

    // copy data from CPU to GPU
    cudaMemcpy((void*)dx, (void*)x, nBytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dy, (void*)y, nBytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)dz, (void*)z, nBytes, cudaMemcpyHostToDevice);

    // kernel configuration
    dim3 blockSize(256);
    dim3 gridSize((N+blockSize.x - 1) / blockSize.x );
    // execute kernel
    add << < gridSize, blockSize >> >(dx, dy, dz, N);

    // copy data from device to host
    cudaMemcpy((void*)z, (void*)dz, nBytes, cudaMemcpyDeviceToHost);

    // check error
    float maxError = 0.0;
    for(int i=0; i<N; i++) {
        maxError = fmax(maxError, fabs(z[i] - 30.0));
    }
    cout  << "Max Error: " << maxError << endl;


    // release host and device memory
    free(x); free(y); free(z);
    cudaFree(dx); cudaFree(dy); cudaFree(dz);

    // time consuming
    gettimeofday( &end, NULL );
    int timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    cout << "time consuming is " << timeuse/1000 << " ms." << endl;
    



    return 0;
}
