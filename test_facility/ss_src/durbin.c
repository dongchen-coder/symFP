
#define N 1024

void durbin(double *y, double* r, double* z) {

	int k, i;
	double alpha, beta, sum;

/*
	y[0] = -r[0];
	beta = 1.0;
	alpha = -r[0];
*/	
	for (k = 1; k < N; k++) {
		beta = (1-alpha*alpha)*beta;
		sum = 0.0;
		for (i=0; i<k; i++) {
			sum += r[k-i-1]*y[i];
		}
		alpha = - (r[k] + sum)/beta;
		
		for (i=0; i<k; i++) {
			z[i] = y[i] + alpha*y[k-i-1];
		}
		for (i=0; i<k; i++) {
			y[i] = z[i];
		}
		y[k] = alpha;
	}	

}


