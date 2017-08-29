#include "../utility/rt.h"

#define NX 1024
#define NY 1024

bool varify(double* A, double* r, double* s, double* p, double* q) {

    double* tempq = (double*)malloc( NX*sizeof(double));
    double* temps = (double*)malloc( NY*sizeof(double));

    for (int i = 0; i< NY; i++) {
        temps[i] = 0.0;
    }

    for (int i = 0; i < NX; i++) {
        tempq[i] = 0.0;
        for (j = 0; j < NY; j++) {
            temps[j] = temps[j] + r[i] * A[i * NY + j];
            tempq[i] = tempq[i] + A[i * NY + j] * p[j];
        }
        if (tempq[i] != q[i]) {
            return false;
        }
    }

    for (int j = 0; j < NY; j++) {
        if (temps[j] != s[j]) {
            return false;
        }
    }

    return true;
}

// void bicg_cpu(int nx, int ny, double* A, double* r, double* s, double* p, double* q) {
//     int i,j;
    
//     for (i = 0; i < NY; i++)
//     {
//         s[i] = 0.0;
//     }
    
//     for (i = 0; i < NX; i++)
//     {
//         q[i] = 0.0;
//         for (j = 0; j < NY; j++)
//         {
//             s[j] = s[j] + r[i] * A[i * NY + j];
//             q[i] = q[i] + A[i * NY + j] * p[j];
//         }
//     }
// }

void bicg_cpu_trace(double* A, double* r, double* s, double* p, double* q, unsigned int nx, int ny) {
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny; j++) {

            s[j] = s[j] + r[i] * A[i * ny + j];
            
            rtTmpAccess(i);                     // load r[i]
            rtTmpAccess(i * ny + j + nx);       // load A[i * ny + j]
            rtTmpAccess(j + (nx * ny) + nx);    // load s[j] 
            rtTmpAccess(j + (ny * ny) + nx);    // store s[j]

            q[i] = q[i] + A[i * ny + j] * p[j];

            rtTmpAccess(i * ny + j + nx);             // load A[i * ny + j]
            rtTmpAccess(j + (ny * nx) + nx + ny);     // load p[j]
            rtTmpAccess(i + (ny * nx) + nx * 2 + ny); // load q[i]
            rtTmpAccess(i + (ny * nx) + nx * 2 + ny); // store q[i]
        }
    }
    return;
}
 


int main() {
    double* A = (double*)malloc( (NX*NY)*sizeof(double));
    double* r = (double*)malloc( NX*sizeof(double));
    double* s = (double*)malloc( NY*sizeof(double));
    double* q = (double*)malloc( NX*sizeof(double));
    double* p = (double*)malloc( NY*sizeof(double));

    for (int i = 0; i < NX; i++) {
        r[i] = i % 256;
    }

    for (int i = 0; i< NY; i++) {
        p[i] = i % 256;
    }

    for (int i = 0; i< NY*NX; i++) {
        A[i] = i % 128;
    }

    bicg_cpu_trace(A, r, s, p ,q, NX, NY);
    dumpRtTmp();

    if (varify(A, r, s, p, q)) {
        cout << "Success" << endl;
    } else {
        cout << "Failed" << endl;
    }

    return 0;
}

