typedef float DATA_TYPE;


	
__kernel void fdtd_kernel1(__global DATA_TYPE *_fict_, __global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int t, int nx, int ny) 
{    
	int j = get_global_id(0);
	int i = get_global_id(1);

	if ((i < nx) && (j < ny))
	{
		int tid = i * ny + j;

		if (i == 0) 
		{
			ey[i * ny + j] = _fict_[t];
		}
		else
		{ 
			ey[i * ny + j] = ey[i * ny + j] - 0.5*(hz[i * ny + j] - hz[(i-1) * ny + j]);
		}
	}
}


__kernel void fdtd_kernel2(__global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int nx, int ny) 
{    
	int j = get_global_id(0);
	int i = get_global_id(1);
	
	if ((i < nx) && (j < ny) && (j > 0))
	{
		ex[i * ny + j] = ex[i * ny + j] - 0.5*(hz[i * ny + j] - hz[i * ny + (j-1)]);
	}
}


__kernel void fdtd_kernel3(__global DATA_TYPE *ex, __global DATA_TYPE *ey, __global DATA_TYPE *hz, int nx, int ny) 
{    
	int j = get_global_id(0);
	int i = get_global_id(1);
	
	if ((i < (nx-1)) && (j < (ny-1)))
	{
		hz[i * ny + j] = hz[i * ny + j] - 0.7*(ex[i * ny + (j+1)] - ex[i * ny + j] + ey[(i + 1) * ny + j] - ey[i * ny + j]);
	}
}
