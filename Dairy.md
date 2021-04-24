# numba GPU 学习内容：

## Intro
使用CUDA时需要"define a thread hierarchy of grid, blocks and threads"。
有三类GPU内存：global device memory， on-chip shared memory and local memory. 使用时需要用户自己决策memory的使用方案来最小化带宽需求和竞争。

## Kernel函数声明：
CPU执行的入口。需要声明线程结构。不能返回值。
CUDA中核函数所定义的运算叫`线程`thread，多个线程组成一个块block，多个块组成网格grid。对应在硬件方面即是，一个thread运行在一个core上，多个thread组成的block运行在SM（streaming multiprocessor）上，多个block组成的grid运行在一个GPU显卡上。因此我们可以通过block和thread号来定位线程。

### block size:
- 软件方面，block size决定了多少个线程共享一块贡献内存；
- 硬件方面， block size要足够大以保证所有执行单元都被占用。
"A thread block size of 16x16 (256 threads), although arbitrary in this case, is a common choice." --from CUDA.
Block中线程的个数一般取32等的倍数，而block的个数则由想要同时并行的线程数目向上取整决定。


### threads positioning
执行kernel时候，kernel函数的代码会被每个线程执行一次。



## 内存管理
Numba会在每次kernel结束后将device memory拷贝到host。为了避免不必要的数据传输，例如对于一些read-only的数组，需要手动控制传输流程：使用cude.to_device来指定拷到device的数据，不需要拷贝回host。

### Example: 矩阵乘法
```
@cuda.jit
def matmul(A, B, C):
    """Perform square matrix multiplication of C = A * B
    """
    i, j = cuda.grid(2)
    if i < C.shape[0] and j < C.shape[1]:
        tmp = 0.
        for k in range(A.shape[1]):
            tmp += A[i, k] * B[k, j]
        C[i, j] = tmp

```
这段代码会重复载入相同的矩阵元素很多次到device memory里。 如果运用blocked algorithm来减少对device memory的访问会提升效率。CUDA能在block里提供一个比较快的shared memory 使得同一个block里的线程可以合作计算一个任务。
```
from numba import cuda, float32

# Controls threads per block and shared memory usage.
# The computation will be done on blocks of TPBxTPB elements.
TPB = 16

@cuda.jit
def fast_matmul(A, B, C):
    # Define an array in the shared memory
    # The size and type of the arrays must be known at compile time
    sA = cuda.shared.array(shape=(TPB, TPB), dtype=float32)
    sB = cuda.shared.array(shape=(TPB, TPB), dtype=float32)

    x, y = cuda.grid(2) # return absolute position of the current thread in the entire grid of blocks...

    tx = cuda.threadIdx.x
    ty = cuda.threadIdx.y
    bpg = cuda.gridDim.x    # blocks per grid

    if x >= C.shape[0] and y >= C.shape[1]:
        # Quit if (x, y) is outside of valid C boundary
        return

    # Each thread computes one element in the result matrix.
    # The dot product is chunked into dot products of TPB-long vectors.
    tmp = 0.
    for i in range(bpg):
        # Preload data into shared memory
        sA[tx, ty] = A[x, ty + i * TPB]
        sB[tx, ty] = B[tx + i * TPB, y]

        # Wait until all threads finish preloading
        cuda.syncthreads()

        # Computes partial product on the shared memory
        for j in range(TPB):
            tmp += sA[tx, j] * sB[j, ty]

        # Wait until all threads finish computing
        cuda.syncthreads()

    C[x, y] = tmp
```




## Device函数：
Device函数只能在GPU内调用，可以返回值。

## Tips：
- CUDA统一的内存管理机制是GPU运行到某处发现数据不再device时，会返回host拷贝数据。当执行完和函数后，又将所有的内存拷贝到host。而一般输入的只读向量是没必要拷贝回host的。
- 流水线优化：一般一边对一些批次的数据进行运算，一般拷贝下一批数据。一个用的是core，一个用的是bus，互不干扰，形成流水线机制。
