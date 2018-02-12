//void atax_cpu(int nx, int ny, DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny), DATA_TYPE POLYBENCH_1D(x,NY,ny), DATA_TYPE POLYBENCH_1D(y,NY,ny), DATA_TYPE POLYBENCH_1D(tmp,NX,nx))

#define NX 1024
#define NY 1024

void atax(int nx, int ny, double* A, double* x, double* y, double* tmp)
{
    int i,j;
    
    for (i= 0; i < NY; i++)
    {
        y[i] = 0;
    }
    
    for (i = 0; i < NX; i++)
    {
      		tmp[i] = 0;
        
      		for (j = 0; j < NY; j++)
            {
                tmp[i] = tmp[i] + A[i * NY + j] * x[j];
            }
        
      		for (j = 0; j < NY; j++)
            {
                y[j] = y[j] + A[i * NY + j] * tmp[i];
            }
    }
    
    return;
}