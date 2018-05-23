
#include "../utility/data_size.h"

#ifdef ORG
    #define NX 1024
    #define NY 1024
#elif defined(TX)
    #define NX 1024
    #define NY 2048
#elif defined(FX)
    #define NX 1024
    #define NY 4096
#elif defined(EX)
    #define NX 1024
    #define NY 8192
#endif

void bicg_cpu(int nx, int ny, double* A, double* r, double* s, double* p, double* q) {
    int i,j;
    
    for (i = 0; i < NY; i++)
    {
        s[i] = 0.0;
    }
    
    for (i = 0; i < NX; i++)
    {
        q[i] = 0.0;
        for (j = 0; j < NY; j++)
        {
            s[j] = s[j] + r[i] * A[i * NY + j];
            q[i] = q[i] + A[i * NY + j] * p[j];
        }
    }
}
