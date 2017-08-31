
int stencil(double* a, double* b) {

	for (int i = 1; i < 1025; i++) {
		for (int j = 1; j < 1025; j++) {
			b[1026*i+j] = a[1026*i+j] 
						+ a[1026*i+j+1] 
						+ a[1026*i+j-1]
						+ a[1026*(i-1)+j]
						+ a[1026*(i+1)+j];
		}
	}

	return 0;
}


/*
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
*/

/*
#define N 1024

void runMvt(int n, double* a, double* x1, double* x2, double* y1, double* y2)
{
    int i, j;

    for (i=0; i< N; i++)
    {
        for (j=0; j < N; j++)
        {
            //x1[i] = x1[i] + a[i][j] * y1[j];
            x1[i] = x1[i] + a[i * N + j] * y1[j];
        }
    }

    for (i=0; i < N; i++)
    {
        for (j=0; j< N; j++)
        {
            //x2[i] = x2[i] + a[j][i] * y2[j];
            x2[i] = x2[i] + a[j * N + i] * y2[j];
        }
    }
}
*/

/*
void test(double* A) {

	for (int i = 0; i < 100; i ++) {
		A[i] = 0;
		for (int j = 0; j < 100; j++) {
			A[i] = 1;
		}
	}

	return;
}
*/
