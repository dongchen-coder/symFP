typedef float DATA_TYPE;


__kernel void doitgen_kernel1(int nr, int nq, int np, __global DATA_TYPE *A, __global DATA_TYPE *C4, __global DATA_TYPE *sum, int r)
{
	int p = get_global_id(0);
	int q = get_global_id(1);

	if ((p < np) && (q < nq))
	{
		sum[r * (nq * np) + q * np + p] = (DATA_TYPE)0.0;
	
		for (int s = 0; s < np; s++)
		{
			sum[r * (nq * np) + q * np + p] = sum[r * (nq * np) + q * np + p] + A[r * (nq * np) + q * np + s] * C4[s * np + p];
		}
	}
}

__kernel void doitgen_kernel2(int nr, int nq, int np, __global DATA_TYPE *A, __global DATA_TYPE *C4, __global DATA_TYPE *sum, int r)
{
	int p = get_global_id(0);
	int q = get_global_id(1);

	if ((p < np) && (q < nq))
	{
		A[r * (nq * np) + q * np + p] = sum[r * (nq * np) + q * np + p];
	}
}
