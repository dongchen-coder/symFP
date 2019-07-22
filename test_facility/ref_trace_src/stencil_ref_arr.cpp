
#include "../utility/rt.h"
#include "../utility/data_size.h"


#ifdef ORG
	#define DIM_SIZE 1024
#elif defined(TX)
	#define DIM_SIZE 1448
#elif defined(FX)
	#define DIM_SIZE 2048
#endif


bool varify(int * b, int * a) {
    for (int i = 1; i < DIM_SIZE+1; i++) {
        for (int j = 1; j < DIM_SIZE+1; j++) {
            if (b[i * (DIM_SIZE+2) + j] != a[i*(DIM_SIZE+2)+j] + a[i*(DIM_SIZE+2)+j + 1] + a[i*(DIM_SIZE+2)+j - 1] + a[(i-1)*(DIM_SIZE+2)+j] + a[(i+1)*(DIM_SIZE+2)+j]) {
                return false;
            }
        }
    }
    return true;
}

void stencil_trace(int *a, int *b, unsigned int dim_size) {
    for (int i = 1; i < dim_size+1; i++) {
        for (int j = 1; j < dim_size+1; j++) {
            b[i* (DIM_SIZE + 2) +j] =  a[i* (DIM_SIZE + 2)+j] + a[i* (DIM_SIZE + 2)+j + 1] + a[i* (DIM_SIZE + 2)+j - 1] + a[(i-1)* (DIM_SIZE + 2) +j] + a[(i+1)* (DIM_SIZE + 2) +j];

            rtTmpAccess(i * (DIM_SIZE + 2) + j, 0, 0);
            rtTmpAccess(i * (DIM_SIZE + 2) + j + 1, 1, 0);
            rtTmpAccess(i * (DIM_SIZE + 2) + j - 1, 2, 0);
            rtTmpAccess( (i-1) * (DIM_SIZE + 2) + j, 3, 0);
            rtTmpAccess( (i+1) * (DIM_SIZE + 2) + j, 4, 0);
            rtTmpAccess(i * (DIM_SIZE + 2) + j + (DIM_SIZE + 2)* (DIM_SIZE + 2), 5, 0);
        }
    }
    return;
}



int main() {
    int* a = (int*)malloc( (DIM_SIZE+2)* (DIM_SIZE+2)*sizeof(int));
    int* b = (int*)malloc( (DIM_SIZE+2)* (DIM_SIZE+2)*sizeof(int));

    for (int i = 0; i < (DIM_SIZE+2) * (DIM_SIZE+2); i++) {
            a[i] = i % 256;
    }

    stencil_trace(a, b, DIM_SIZE);
    
    dumpRI();
    /*
    if (varify(b, a)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }
     */

    return 0;
}
