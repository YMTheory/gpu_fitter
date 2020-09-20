#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

using namespace std;

int main()
{
    struct timeval start, end;
    gettimeofday( &start, NULL );
    float *A, *B, *C;
    int n = 1024 * 1024;
    int size = n * sizeof(float);
    A = (float*)malloc(size);
    B = (float*)malloc(size);
    C = (float*)malloc(size);

    for(int i=0; i<n; i++)
    {
        A[i] = 90.0;
        B[i] = 10.0;
    }
    
    for(int i=0; i<n; i++) {
        C[i] = A[i] + B[i];
    }

    float max_error = 0.0;
    for(int i=0; i<n; i++) {
        max_error += fabs(100-C[i]);
    }

    cout << "max error is : " << max_error << endl;
    gettimeofday( &end, NULL );
    int timeuse = 100000 * (end.tv_sec - start.tv_sec ) + end.tv_usec - start.tv_usec;
    cout << "total time consumed is : " << timeuse / 10000 << "ms" << endl;
    
    return 0;

}

