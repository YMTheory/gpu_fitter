from numba import cuda
import numpy as np

@cuda.jit
def increment_by_one(an_array):
    """
    Increment all array elements by one
    """
    # thread positioning
    tx = cuda.threadIdx.x
    ty = cuda.blockIdx.x
    bw = cuda.blockDim.x
    pos = tx + ty * bw
    if pos < an_array.size: 
        an_array[pos] += 1


def main():
    arr = np.arange(0, 10, 1)
    threadsperblock = 32
    blockspergrid = (arr.size + (threadsperblock - 1))
    increment_by_one[blockspergrid, threadsperblock](arr)
    cuda.synchronize()
    print(arr)



if __name__ == "__main__":
    main()