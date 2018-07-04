typedef float DATA_TYPE;



__kernel void adi_kernel1(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int n)
{
	int i1 = get_global_id(0);
	int i2;	

	if ((i1 < n))
	{
		for (i2 = 1; i2 < n; i2++)
		{
			X[i1*n + i2] = X[i1*n + i2] - X[i1*n + (i2-1)] * A[i1*n + i2] / B[i1*n + (i2-1)];
			B[i1*n + i2] = B[i1*n + i2] - A[i1*n + i2] * A[i1*n + i2] / B[i1*n + (i2-1)];
		}
	}
}

__kernel void adi_kernel2(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int n)
{
	int i1 = get_global_id(0);
	
	if ((i1 < n))
	{
		X[i1*n + (n-1)] = X[i1*n + (n-1)] / B[i1*n + (n-1)];
	}
}
	
__kernel void adi_kernel3(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int n)
{
	int i1 = get_global_id(0);
	int i2;	

	if ((i1 < n))
	{
		for (i2 = 0; i2 < n-2; i2++)
		{
			X[i1*n + (n-i2-2)] = (X[i1*n + (n-2-i2)] - X[i1*n + (n-2-i2-1)] * A[i1*n + (n-i2-3)]) / B[i1*n + (n-3-i2)];
		}
	}
}



__kernel void adi_kernel4(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int i1, int n)
{
	int i2 = get_global_id(0);
	
	if ((i2 < n))
	{
		X[i1*n + i2] = X[i1*n + i2] - X[(i1-1)*n + i2] * A[i1*n + i2] / B[(i1-1)*n + i2];
		B[i1*n + i2] = B[i1*n + i2] - A[i1*n + i2] * A[i1*n + i2] / B[(i1-1)*n + i2];
	}
}

__kernel void adi_kernel5(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int n)
{
	int i2 = get_global_id(0);
	
	if ((i2 < n))
	{
		X[(n-1)*n + i2] = X[(n-1)*n + i2] / B[(n-1)*n + i2];
	}
}

__kernel void adi_kernel6(__global DATA_TYPE* A, __global DATA_TYPE* B, __global DATA_TYPE* X, int i1, int n)
{
	int i2 = get_global_id(0);
	
	if ((i2 < n))
	{
	     X[(n-2-i1)*n + i2] = (X[(n-2-i1)*n + i2] - X[(n-i1-3)*n + i2] * A[(n-3-i1)*n + i2]) / B[(n-2-i1)*n + i2];
	}
}
