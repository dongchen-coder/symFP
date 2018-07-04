typedef float DATA_TYPE;



__kernel void mvt_kernel1(__global DATA_TYPE *a, __global DATA_TYPE *x1, __global DATA_TYPE *y1, int n) 
{    
	int i = get_global_id(0);

	if (i < n)
	{
		int j;	
		for (j=0; j < n; j++)
		{
			x1[i] += a[i * n + j] * y1[j];
		}
	}
}

__kernel void mvt_kernel2(__global DATA_TYPE *a, __global DATA_TYPE *x2, __global DATA_TYPE *y2, int n) 
{    
	int i = get_global_id(0);

	if (i < n)
	{
		int j;	
		for (j=0; j < n; j++)
		{
			x2[i] += a[j * n + i] * y2[j];	
		}
	}
}
