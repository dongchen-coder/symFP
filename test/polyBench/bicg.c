//void bicg_cpu(int nx, int ny, DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny), DATA_TYPE POLYBENCH_1D(r,NX,nx), DATA_TYPE POLYBENCH_1D(s,NY,ny), DATA_TYPE POLYBENCH_1D(p,NY,ny), DATA_TYPE POLYBENCH_1D(q,NX,nx))

#define NX 1024
#define NY 1024

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
