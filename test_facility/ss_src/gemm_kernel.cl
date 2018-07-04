typedef float DATA_TYPE;

__kernel void gemm(__global DATA_TYPE *a, __global DATA_TYPE *b, __global DATA_TYPE *c, DATA_TYPE alpha, DATA_TYPE beta, int ni, int nj, int nk) 
{
    	int j = get_global_id(0);
	int i = get_global_id(1);

	if ((i < ni) && (j < nj))
	{	
		c[i * nj + j] *= beta;
		int k;
		for(k=0; k < nk; k++)
		{
			c[i * nj + j] += alpha * a[i * nk + k] * b[k * nj +j];
		}
	}
}
