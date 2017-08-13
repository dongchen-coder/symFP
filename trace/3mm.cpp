#include "./utility/rt.h"

#define DIM_SIZE 1024

bool varify(int *G, int *A, int *B, int *C, int *D, int *E, int *F) {
    int i, j, k, temp;
    
    /* E := A*B */
    for (i = 0; i < DIM_SIZE; i++)
    {
        for (j = 0; j < DIM_SIZE; j++)
        {
            temp = 0;
            for (k = 0; k < DIM_SIZE; ++k)
            {
                temp += A[i * DIM_SIZE + k] * B[k * DIM_SIZE + j];
            }
            if (E[i * DIM_SIZE + j] != temp) {
                return false;
            }
        }
    }
    
    /* F := C*D */
    for (i = 0; i < DIM_SIZE; i++)
    {
        for (j = 0; j < DIM_SIZE; j++)
        {
            temp = 0;
            for (k = 0; k < DIM_SIZE; ++k)
            {
                temp += C[i * DIM_SIZE + k] * D[k * DIM_SIZE + j];
            }
            if (F[i * DIM_SIZE + j] != temp) {
                return false;
            }
        }
    }
    
    /* G := E*F */
    for (i = 0; i < DIM_SIZE; i++)
    {
        for (j = 0; j < DIM_SIZE; j++)
        {
            temp = 0;
            for (k = 0; k < DIM_SIZE; ++k)
            {
                temp += E[i * DIM_SIZE + k] * F[k * DIM_SIZE + j];
            }
            if (G[i * DIM_SIZE + j] != temp) {
                return false;
            }
        }
    }

    return true;
}

// void mm3_cpu(int ni, int nj, int nk, int nl, int nm,
//              double * E, double* A, double* B, double* F, double* C, double* D, double* G)
// {
//     int i, j, k;
    
//     /* E := A*B */
//     for (i = 0; i < NI; i++)
//     {
//         for (j = 0; j < NJ; j++)
//         {
//             E[i * NJ + j] = 0;
//             for (k = 0; k < NK; ++k)
//             {
//                 E[i * NJ + j] += A[i * NK + k] * B[k * NJ + j];
//             }
//         }
//     }
    
//     /* F := C*D */
//     for (i = 0; i < NJ; i++)
//     {
//         for (j = 0; j < NL; j++)
//         {
//             F[i * NL + j] = 0;
//             for (k = 0; k < NM; ++k)
//             {
//                 F[i * NL + j] += C[i * NM + k] * D[k * NL + j];
//             }
//         }
//     }
    
//     /* G := E*F */
//     for (i = 0; i < NI; i++)
//     {
//         for (j = 0; j < NL; j++)
//         {
//             G[i * NL + j] = 0;
//             for (k = 0; k < NJ; ++k)
//             {
//                 G[i * NL + j] += E[i * NJ + k] * F[k * NL + j];
//             }
//         }
//     }
// }


void mm3_cpu_trace(int *A, int *B, int *C, int *D, int *E, int *F, int *G, unsigned int dim_size) {
    /* E := A*B */
    for (int i = 0; i < dim_size; i++)
    {
        for (int j = 0; j < dim_size; j++)
        {
            E[i * dim_size + j] = 0;
            for (int k = 0; k < dim_size; ++k)
            {
                E[i * dim_size + j] = E[i * dim_size + j] + A[i * dim_size + k] * B[k * dim_size + j];
                rtTmpAccess(i * dim_size + k);                               // load A[i * NK + k]
                rtTmpAccess(k * dim_size + j + (dim_size * dim_size));       // load B[k * NJ + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 2);   // load E[i * NJ + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 2);   // store E[i * NJ + j]
            }
        }
    }
    
    /* F := C*D */
    for (int i = 0; i < dim_size; i++)
    {
        for (int j = 0; j < dim_size; j++)
        {
            F[i * dim_size + j] = 0;
            for (int k = 0; k < dim_size; ++k)
            {
                F[i * dim_size + j] = F[i * dim_size + j] + C[i * dim_size + k] * D[k * dim_size + j];
                rtTmpAccess(i * dim_size + k + (dim_size * dim_size) * 3);      // load C[i * dim_size + k]
                rtTmpAccess(k * dim_size + j + (dim_size * dim_size) * 4);      // load D[k * dim_size + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 5);      // load F[i * dim_size + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 5);      // stroe F[i * dim_size + j]
            }
        }
    }
    
    /* G := E*F */
    for (int i = 0; i < dim_size; i++)
    {
        for (int j = 0; j < dim_size; j++)
        {
            G[i * dim_size + j] = 0;
            for (int k = 0; k < dim_size; ++k)
            {
                G[i * dim_size + j] = G[i * dim_size + j] + E[i * dim_size + k] * F[k * dim_size + j];
                rtTmpAccess(i * dim_size + k * (dim_size * dim_size) * 2);      // load E[i * dim_size + k]
                rtTmpAccess(k * dim_size + j * (dim_size * dim_size) * 5);      // load F[k * dim_size + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 6);      // load G[i * dim_size + j]
                rtTmpAccess(i * dim_size + j + (dim_size * dim_size) * 6);      // store G[i * dim_size + j]
            }
        }
    }
    return;
}



int main() {
    int* a = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* b = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* c = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* d = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* e = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* f = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );
    int* g = (int*)malloc( (DIM_SIZE)* (DIM_SIZE)*sizeof(int) );


    for (int i = 0; i < (DIM_SIZE * DIM_SIZE); i++) {
            a[i] = i % 256;
            b[i] = i % 128;
    }

    mm3_cpu_trace(a, b, c, d, e, f, g, DIM_SIZE);
    dumpRtTmp();

    if (varify(g, a, b, c, d, e, f)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
}

