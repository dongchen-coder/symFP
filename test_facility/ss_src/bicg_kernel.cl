typedef float DATA_TYPE;

__kernel void bicgKernel1(__global DATA_TYPE *A, __global DATA_TYPE *p, __global DATA_TYPE *q, int nx, int ny) 
{
    	int i = get_global_id(0);
	
	if (i < nx)
	{
		q[i] = 0.0;

		int j;
		for(j=0; j < ny; j++)
		{
			q[i] += A[i * ny + j] * p[j];
		}
	}
	
}

__kernel void bicgKernel2(__global DATA_TYPE *A, __global DATA_TYPE *r, __global DATA_TYPE *s, int nx, int ny) 
{
	int j = get_global_id(0);
	
	if (j < ny)
	{
		s[j] = 0.0;

		int i;
		for(i = 0; i < nx; i++)
		{
			s[j] += A[i * ny + j] * r[i];
		}
	}
	
}
