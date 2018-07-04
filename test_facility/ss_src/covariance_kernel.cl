typedef float DATA_TYPE;


__kernel void mean_kernel(__global DATA_TYPE *mean, __global DATA_TYPE *data, DATA_TYPE float_n, int m, int n) 
{    	 
	int j = get_global_id(0);
	
	if (j < m)
	{
		mean[j] = 0.0;

		int i;
		for(i = 0; i < n; i++)
		{
			mean[j] += data[i * m + j];
		}
		mean[j] /= (DATA_TYPE)float_n;
	}
}

__kernel void reduce_kernel(__global DATA_TYPE *mean, __global DATA_TYPE *data, int m, int n) 
{
	int j = get_global_id(0);    
	int i = get_global_id(1);

	if ((i < n) && (j < m))
	{
		data[i * m + j] -= mean[j];	
	}
}


__kernel void covar_kernel(__global DATA_TYPE *symmat, __global DATA_TYPE *data, int m, int n) 
{
	int j1 = get_global_id(0);
	int i, j2;

	if (j1 < m)
	{
		for (j2 = j1; j2 < m; j2++)
		{		
	      		symmat[j1*m + j2] = 0.0;
			for(i = 0; i < n; i++)
			{
				symmat[j1 * m + j2] += data[i * m + j1] * data[i * m + j2];
			}
			symmat[j2 * m + j1] = symmat[j1 * m + j2];
		}
	}
}
