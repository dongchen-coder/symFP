typedef float DATA_TYPE;



__kernel void lu_kernel1(__global DATA_TYPE *A, int k, int n)
{
	int j = get_global_id(0) + (k + 1);
	
	if ((j < n))
	{
		A[k*n + j] = A[k*n + j] / A[k*n + k];
	}
}

__kernel void lu_kernel2(__global DATA_TYPE *A, int k, int n)
{
	int j = get_global_id(0) + (k + 1);
	int i = get_global_id(1) + (k + 1);
	
	if ((i < n) && (j < n))
	{
		A[i*n + j] = A[i*n + j] - A[i*n + k] * A[k*n + j];
	}
}
