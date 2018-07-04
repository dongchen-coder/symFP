typedef float DATA_TYPE;
	
__kernel void gemver_kernel1(__global DATA_TYPE *A, __global DATA_TYPE *V1, __global DATA_TYPE *V2, __global DATA_TYPE *U1, __global DATA_TYPE *U2, int n) 
{    
	int j = get_global_id(0);
	int i = get_global_id(1);
	
	if ((i < n) && (j < n))
	{
		A[i*n + j] += U1[i] * V1[j] + U2[i] * V2[j];
	}
}


__kernel void gemver_kernel2(__global DATA_TYPE *A, __global DATA_TYPE *X, __global DATA_TYPE *Y, __global DATA_TYPE *Z, DATA_TYPE beta, int n) 
{    
	int i = get_global_id(0);

	if (i < n)
	{
		int j;
		for(j = 0; j < n; j++) 
		{
			X[i] += beta * A[j * n + i] * Y[j];
		}
		X[i] += Z[i];
	}
}


__kernel void gemver_kernel3(__global DATA_TYPE *A, __global DATA_TYPE *X, __global DATA_TYPE *w, DATA_TYPE alpha, int n) 
{    
	int i = get_global_id(0);
	
	if (i < n)
	{
		int j;
		for(j = 0; j < n; j++)
		{ 
			w[i] += alpha * A[i*n + j] * X[j];
		}
	}
}
