
int stencil(double* a, double* b) {

	for (int i = 1; i < 1025; i++) {
		for (int j = 1; j < 1025; j++) {
			b[1026*i+j] = a[1026*i+j] 
						+ a[1026*i+j+1] 
						+ a[1026*i+j-1]
						+ a[1026*(i-1)+j]
						+ a[1026*(i+1)+j];
		}
	}

	return 0;
}

