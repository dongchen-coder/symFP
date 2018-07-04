typedef float DATA_TYPE;

__kernel void atax_kernel1(__global DATA_TYPE *A, __global DATA_TYPE *x, __global DATA_TYPE *tmp, int nx, int ny) {
    
	int i = get_global_id(0);

	if (i < nx)
	{
		int j;
		for(j=0; j < ny; j++)
		{
			tmp[i] += A[i * ny + j] * x[j];
		}
	}
}

__kernel void atax_kernel2(__global DATA_TYPE *A, __global DATA_TYPE *y, __global DATA_TYPE *tmp, int nx, int ny) {
    
	int j = get_global_id(0);

	if (j < ny)
	{
		int i;
		for(i=0; i < nx; i++)
		{
			y[j] += A[i * ny + j] * tmp[i];
		}
	}
}
