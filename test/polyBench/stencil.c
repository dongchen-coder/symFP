


/*
int stencil(double* a, double* b) {

	for (int i = 0; i < 1025; i++) {
		for (int j = 0; j < 1025; j++) {
			b[1026*i+j] = a[1026*i+j] 
						+ a[1026*i+j+1] 
						+ a[1026*i+j-1]
						+ a[1026*(i+1)+j]
						+ a[1026*(i-1)+j];
		}
	}

	return 0;
}
*/

#define NX 1024
#define NY 1024

void atax_cpu(int nx, int ny, double* A, double* x, double* y, double* tmp)
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

