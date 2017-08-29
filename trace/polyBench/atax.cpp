#include "../utility/rt.h"

#define NX 1024
#define NY 1024

bool varify(double* A, double* x, double* y, double* tmp) {

    // init temp array
    double* temptmp = (double*)malloc( NX*sizeof(double) );
    double* tempy = (double*)malloc( NY*sizeof(double) );

    for (i = 0; i< NY; i++) {
        tempy[i] = 0;
    }

    for (i = 0; i < NX; i++) {

        temptmp[i] = 0;
    
        for (j = 0; j < NY; j++) {
            temptmp[i] = temptmp[i] + A[i * NY + j] * x[j];
        }
        if (temptmp[i] != tmp[i]) {
            return false;
        }
    
        for (j = 0; j < NY; j++) {
            tempy[j] = tempy[j] + A[i * NY + j] * tmp[i];
        }
    }

    for (i = 0; i< NY; i++) {
        if (tempy[i] != y[i]) {
            return false;
        }
    }
    return true;
}

// #define NX 1024
// #define NY 1024

// void atax_cpu(int nx, int ny, double* A, double* x, double* y, double* tmp)
// {
//     int i,j;
    
//     for (i= 0; i < NY; i++)
//     {
//         y[i] = 0;
//     }
    
//     for (i = 0; i < NX; i++)
//     {
//             tmp[i] = 0;
        
//             for (j = 0; j < NY; j++)
//             {
//                 tmp[i] = tmp[i] + A[i * NY + j] * x[j];
//             }
        
//             for (j = 0; j < NY; j++)
//             {
//                 y[j] = y[j] + A[i * NY + j] * tmp[i];
//             }
//     }
    
//     return;
// }

void atax_cpu_trace(double* A, double* x, double* y, double* tmp, unsigned int dim_size_x, unsigned int dim_size_y) {
    int i,j;
    
    for (i= 0; i < dim_size_y; i++)
    {
        y[i] = 0;
    }
    
    for (i = 0; i < dim_size_x; i++)
    {
            tmp[i] = 0;
        
            for (j = 0; j < dim_size_y; j++)
            {
                tmp[i] = tmp[i] + A[i * dim_size_y + j] * x[j];
                rtTmpAccess(i * dim_size_y + j);                // load A[i * dim_size_y + j]
                rtTmpAccess(j + dim_size_x * dim_size_y);       // load x[j]
                rtTmpAccess(i + dim_size_x * dim_size_y * 2);   // load tmp[i]
                rtTmpAccess(i + dim_size_x * dim_size_y * 2);   // store tmp[i]    
            }
        
            for (j = 0; j < dim_size_y; j++)
            {
                y[j] = y[j] + A[i * dim_size_y + j] * tmp[i];
                rtTmpAccess(i * dim_size_y + j);                // load A[i * dim_size_y + j]
                rtTmpAccess(i + dim_size_x * dim_size_y * 2);   // load tmp[i]
                rtTmpAccess(i + dim_size_x * dim_size_y * 3);   // load y[j]
                rtTmpAccess(i + dim_size_x * dim_size_y * 3);   // store y[j]

            }
    }
    return;
}



int main() {

    double* A = (double*)malloc( (NX * NY) * sizeof(double) );
    double* x = (double*)malloc( NX * sizeof(double) );
    double* y = (double*)malloc( NY * sizeof(double) );
    double* tmp = (double*)malloc( NX * sizeof(double) );

    for (int i = 0; i < (NX * NY); i++) {
        A[i] = i % 256;
    }

    for (int i = 0; i < NX; i++) {
        x[i] = i % 256;
    }

    atax_cpu_trace(A, x, y, tmp, NX, NY);
    dumpRtTmp();

    if (varify(A, x, y, tmp)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
}

