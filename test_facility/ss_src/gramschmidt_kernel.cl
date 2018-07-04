typedef float DATA_TYPE;


__kernel void gramschmidt_kernel1(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int ni, int nj)
{
	int tid = get_global_id(0);
	
	if (tid == 0)
	{
		DATA_TYPE nrm = 0.0;
		int i;
		for (i = 0; i < ni; i++)
		{
			nrm += a[i * nj + k] * a[i * nj + k];
		}
      		r[k * nj + k] = sqrt(nrm);
	}
}


__kernel void gramschmidt_kernel2(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int ni, int nj)
{
	int i = get_global_id(0);

        if (i < ni)
	{	
		q[i * nj + k] = a[i * nj + k] / r[k * nj + k];
	}
}


__kernel void gramschmidt_kernel3(__global DATA_TYPE *a, __global DATA_TYPE *r, __global DATA_TYPE *q, int k, int ni, int nj)
{
	int j = get_global_id(0) + (k+1);

	if ((j < nj))
	{
		r[k*nj + j] = 0.0;

		int i;
		for (i = 0; i < ni; i++)
		{
			r[k*nj + j] += q[i*nj + k] * a[i*nj + j];
		}
		
		for (i = 0; i < ni; i++)
		{
			a[i*nj + j] -= q[i*nj + k] * r[k*nj + j];
		}
	}
}
