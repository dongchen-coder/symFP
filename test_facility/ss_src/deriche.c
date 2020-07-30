

#ifndef DEBUG
# if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(MEDIUM_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
#  define LARGE_DATASET
# endif
#ifdef MINI_DATASET
    #define W 32
    #define H 32
#endif
#ifdef SMALL_DATASET
    #define W 1024
    #define H 1024
#endif 
#ifdef MEDIUM_DATASET
    #define W 4096
    #define H 4096
#endif
#ifdef LARGE_DATASET
    #define W 8192
    #define H 8192
#endif
#ifdef EXTRALARGE_DATASET
    #define W 100000
    #define H 100000
#endif

#else
    #define W 4
    #define H 4
#endif

#define EXP_FUN(x) exp(x)
#define POW_FUN(x,y) pow(x,y)

void deriche(double* y1, double* imgIn, double* y2, double* imgOut, double alpha) {

	double k, a1, a2, a3, a4, b1, b2, c1, c2;

	double a5, a6, a7, a8;

	double ym1, ym2, xm1, yp1, yp2, xp1, xp2;

	double tm1;

	double tp1, tp2;

	int i, j;

	k = (1.0 - EXP_FUN(-alpha))*(1.0-EXP_FUN(-alpha))/(1.0+2.0*alpha*EXP_FUN(-alpha)-EXP_FUN(2.0*alpha));
	a1 = a5 = k;
	a2 = a6 = k*EXP_FUN(-alpha)*(alpha-1.0);
	a3 = a7 = k*EXP_FUN(-alpha)*(alpha+1.0);
	a4 = a8 = -k*EXP_FUN(-2.0*alpha);
	b1 =  POW_FUN(2.0,-alpha);
	b2 = -EXP_FUN(-2.0*alpha);
	c1 = c2 = 1;

	for (i=0; i < W; i++) {
		ym1 = 0.0;
		ym2 = 0.0;
		xm1 = 0.0;
		for (j=0; j < H; j++) {
			y1[i * H + j] = a1*imgIn[i * H + j] + a2*xm1 + b1*ym1 + b2*ym2;
			xm1 = imgIn[i * H + j];
			ym2 = ym1;
			ym1 = y1[i * H + j];
		}
	}

	for (i=0; i < W; i++) {
		yp1 = 0.0;
		yp2 = 0.0;
		xp1 = 0.0;
		xp2 = 0.0;
		for (j = 0; j < H; j++) {
			y2[i * H + H-1-j] = a3*xp1 + a4*xp2 + b1*yp1 + b2*yp2;
			xp2 = xp1;
			xp1 = imgIn[i * H + H-1-j];
			yp2 = yp1;
			yp1 = y2[i * H + H-1-j];
		}
	}

	for (i=0; i < W; i++)
		for (j=0; j < H; j++) {
			imgOut[i * H + j] = c1 * (y1[i * H + j] + y2[i * H + j]);
		}

	for (j=0; j < H; j++) {
		tm1 = 0.0;
		ym1 = 0.0;
		ym2 = 0.0;
		for (i=0; i < W; i++) {
			y1[i * H + j] = a5*imgOut[i * H + j] + a6*tm1 + b1*ym1 + b2*ym2;
			tm1 = imgOut[i * H + j];
			ym2 = ym1;
			ym1 = y1 [i * H + j];
		}
	}

	for (j=0; j < H; j++) {
		tp1 = 0.0;
		tp2 = 0.0;
		yp1 = 0.0;
		yp2 = 0.0;
		for (i = 0; i < W; i++) {
			y2[(W-1-i) * H + j] = a7*tp1 + a8*tp2 + b1*yp1 + b2*yp2;
			tp2 = tp1;
			tp1 = imgOut[(W-1-i) * H + j];
			yp2 = yp1;
			yp1 = y2[(W-1-i) * H + j];
		}
	}

	for (i=0; i < W; i++)
		for (j=0; j < H; j++)
			imgOut[i * H + j] = c2*(y1[i * H + j] + y2[i * H + j]);

}

