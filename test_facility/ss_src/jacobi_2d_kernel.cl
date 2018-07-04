typedef float DATA_TYPE;


__kernel void runJacobi2D_kernel1(__global DATA_TYPE* A, __global DATA_TYPE* B, int n)
{
	int i = get_global_id(1);
	int j = get_global_id(0);

	if ((i >= 1) && (i < (n-1)) && (j >= 1) && (j < (n-1)))
	{
		B[i*n + j] = 0.2f * (A[i*n + j] + A[i*n + (j-1)] + A[i*n + (1 + j)] + A[(1 + i)*n + j] 
				+ A[(i-1)*n + j]);	
	}
}

__kernel void runJacobi2D_kernel2(__global DATA_TYPE* A, __global DATA_TYPE* B, int n)
{
	int i = get_global_id(1);
	int j = get_global_id(0);
	
	if ((i >= 1) && (i < (n-1)) && (j >= 1) && (j < (n-1)))
	{
		A[i*n + j] = B[i*n + j];
	}
}
