typedef float DATA_TYPE;
	
__kernel void gesummv_kernel(__global DATA_TYPE *a, __global DATA_TYPE *b, __global DATA_TYPE *x, __global DATA_TYPE *y, __global DATA_TYPE *tmp, DATA_TYPE alpha, DATA_TYPE beta, int n) 
{    
	int i = get_global_id(0);

	if (i < n)
	{
		int j;
		for(j = 0; j < n; j++)
		{	
			tmp[i] += a[i * n + j] * x[j];
			y[i] += b[i * n + j] * x[j];
		}
		y[i] = alpha * tmp[i] + beta * y[i];
	}
}
