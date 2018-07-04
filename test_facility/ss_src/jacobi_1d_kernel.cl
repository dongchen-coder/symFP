typedef float DATA_TYPE;

__kernel void runJacobi1D_kernel1(__global DATA_TYPE* A, __global DATA_TYPE* B, int n)
{
	int i = get_global_id(0);
	if ((i >= 1) && (i < (n-1)))
	{
		B[i] = 0.33333f * (A[i-1] + A[i] + A[i + 1]);
	}
}

__kernel void runJacobi1D_kernel2(__global DATA_TYPE* A, __global DATA_TYPE* B, int n)
{
	int j = get_global_id(0);
	
	if ((j >= 1) && (j < (n-1)))
	{
		A[j] = B[j];
	}
}
