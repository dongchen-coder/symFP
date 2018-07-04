typedef float DATA_TYPE;


/*  C := alpha*A*A' + beta*C */
__kernel void syrk_kernel(__global DATA_TYPE *a, __global DATA_TYPE *c, DATA_TYPE alpha, DATA_TYPE beta, int ni, int nj) 
{
   	int j = get_global_id(0);
	int i = get_global_id(1);
	
	if ((i < nj) && (j < nj))
	{
		c[i * nj + j] *= beta;
		int k;		
		for(k=0; k< ni; k++)
		{
			c[i * nj + j] += alpha * a[i * ni + k] * a[j * ni + k];
		}
	}
}
