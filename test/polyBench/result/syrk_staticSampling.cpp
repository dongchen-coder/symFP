
 /* Start to analysis array index
Array index info
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)
A.addr ((i * 256) + k)
A.addr ((j * 256) + k)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
double %alpha
double %beta
double* %A
double* %C

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access   %6 = load double, double* %arrayidx, align 8
------array access   store double %mul4, double* %arrayidx, align 8
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------k
------Loop Bound: (0, 256)
--------array access   %16 = load double, double* %arrayidx20, align 8
--------array access   %20 = load double, double* %arrayidx25, align 8
--------array access   %24 = load double, double* %arrayidx30, align 8
--------array access   store double %add31, double* %arrayidx30, align 8

Finish analysis loops */ 
 // Start to generating Static Sampling Code
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
using namespace std;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( int rt) {
    if (RT.find(rt) == RT.end()) { 
        RT[rt] = 1;
    } else {
        RT[rt] += 1;
    }
    return;
}
void RTtoMR_AET() {
    std::map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        P[it->first] = accumulate_num_RT / total_num_RT;
        accumulate_num_RT += it->second;
    }
    P[0] = 1;
    double sum_P = 0;
    uint64_t t = 0;
    uint64_t prev_t = 0;
    for (uint64_t c = 0; c <= max_RT; c++) {
        while (sum_P < c && t <= max_RT) {
            if (P.find(t) != P.end()) {
                sum_P += P[t];
                prev_t = t;
            } else {
                sum_P += P[prev_t];
            }
            t++;
        }
        MR[c] = P[prev_t];
    }
    return;
}
int calAddrC_addr0( int i, int j) {
    int result = ((i * 256) + j);
    return result;
}
int calAddrC_addr1( int i, int j) {
    int result = ((i * 256) + j);
    return result;
}
int calAddrA_addr0( int i, int j, int k) {
    int result = ((i * 256) + k);
    return result;
}
int calAddrA_addr1( int i, int j, int k) {
    int result = ((j * 256) + k);
    return result;
}
int calAddrC_addr2( int i, int j, int k) {
    int result = ((i * 256) + j);
    return result;
}
int calAddrC_addr3( int i, int j, int k) {
    int result = ((i * 256) + j);
    return result;
}
int rtCalA_addr0_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 0 - 0;
}
int rtCalA_addr0_1(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 1 - 0;
}
int rtCalA_addr1_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 0 - 1;
}
int rtCalA_addr1_1(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 1 - 1;
}
int rtCalC_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 512 + (jreuse - j) * 2 + 0 - 0;
}
int rtCalC_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 512 + (jreuse - j) * 2 + 1 - 0;
}
int rtCalC_addr0_2(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (256 - i) * 512 + (ireuse - 0) * 262144 + (256 - j) * 2 + (jreuse - 0) * 1024 + (kreuse - 0) * 4 + 2 - 0;
}
int rtCalC_addr0_3(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (256 - i) * 512 + (ireuse - 0) * 262144 + (256 - j) * 2 + (jreuse - 0) * 1024 + (kreuse - 0) * 4 + 3 - 0;
}
int rtCalC_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 512 + (jreuse - j) * 2 + 0 - 1;
}
int rtCalC_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 512 + (jreuse - j) * 2 + 1 - 1;
}
int rtCalC_addr1_2(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (256 - i) * 512 + (ireuse - 0) * 262144 + (256 - j) * 2 + (jreuse - 0) * 1024 + (kreuse - 0) * 4 + 2 - 1;
}
int rtCalC_addr1_3(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (256 - i) * 512 + (ireuse - 0) * 262144 + (256 - j) * 2 + (jreuse - 0) * 1024 + (kreuse - 0) * 4 + 3 - 1;
}
int rtCalC_addr2_0(int i, int j, int k, int ireuse, int jreuse) {
    return (256 - i) * 262144 + (ireuse - 0) * 512 + (256 - j) * 1024 + (jreuse - 0) * 2 + (256 - k) * 4 + 0 - 2;
}
int rtCalC_addr2_1(int i, int j, int k, int ireuse, int jreuse) {
    return (256 - i) * 262144 + (ireuse - 0) * 512 + (256 - j) * 1024 + (jreuse - 0) * 2 + (256 - k) * 4 + 1 - 2;
}
int rtCalC_addr2_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 2 - 2;
}
int rtCalC_addr2_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 3 - 2;
}
int rtCalC_addr3_0(int i, int j, int k, int ireuse, int jreuse) {
    return (256 - i) * 262144 + (ireuse - 0) * 512 + (256 - j) * 1024 + (jreuse - 0) * 2 + (256 - k) * 4 + 0 - 3;
}
int rtCalC_addr3_1(int i, int j, int k, int ireuse, int jreuse) {
    return (256 - i) * 262144 + (ireuse - 0) * 512 + (256 - j) * 1024 + (jreuse - 0) * 2 + (256 - k) * 4 + 1 - 3;
}
int rtCalC_addr3_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 2 - 3;
}
int rtCalC_addr3_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262144 + (jreuse - j) * 1024 + (kreuse - k) * 4 + 3 - 3;
}
bool checkIntervenA_addr0_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrA_addr1(iInterven, jInterven, kInterven) == calAddrA_addr0(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenA_addr0_1(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrA_addr0(iInterven, jInterven, kInterven) == calAddrA_addr0(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenA_addr1_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrA_addr1(iInterven, jInterven, kInterven) == calAddrA_addr1(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenA_addr1_1(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrA_addr0(iInterven, jInterven, kInterven) == calAddrA_addr1(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr0_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenC_addr0_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenC_addr0_2(int i, int j, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == 0 && jInterven == 0) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr3(iInterven, jInterven, kInterven) == calAddrC_addr0(i, j)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr0_3(int i, int j, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == 0 && jInterven == 0) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr2(iInterven, jInterven, kInterven) == calAddrC_addr0(i, j)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr1_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenC_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenC_addr1_2(int i, int j, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == 0 && jInterven == 0) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr3(iInterven, jInterven, kInterven) == calAddrC_addr1(i, j)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr1_3(int i, int j, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 256; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 256) {
            jIntervenUB = 256;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == 0 && jInterven == 0) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr2(iInterven, jInterven, kInterven) == calAddrC_addr1(i, j)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr2_0(int i, int j, int k, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenC_addr2_1(int i, int j, int k, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenC_addr2_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr3(iInterven, jInterven, kInterven) == calAddrC_addr2(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr2_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr2(iInterven, jInterven, kInterven) == calAddrC_addr2(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr3_0(int i, int j, int k, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenC_addr3_1(int i, int j, int k, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenC_addr3_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse - 1 ;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr3(iInterven, jInterven, kInterven) == calAddrC_addr3(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
bool checkIntervenC_addr3_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 256- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse) {
                kIntervenUB = kreuse;
            } else {
                kIntervenUB = 256- 1;
            }
            for(; kInterven <= kIntervenUB; kInterven++) {
                if( calAddrC_addr2(iInterven, jInterven, kInterven) == calAddrC_addr3(i, j, k)) {
                    return true;
                }
            }
        }
    }
    return false;
}
void pairA_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrA_addr0(i, j, k) == calAddrA_addr0(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenA_addr0_0(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalA_addr0_0(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairA_addr0_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrA_addr0(i, j, k) == calAddrA_addr1(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenA_addr0_1(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalA_addr0_1(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairA_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrA_addr1(i, j, k) == calAddrA_addr0(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenA_addr1_0(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalA_addr1_0(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairA_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrA_addr1(i, j, k) == calAddrA_addr1(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenA_addr1_1(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalA_addr1_1(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairC_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                    if ( calAddrC_addr0(i, j) == calAddrC_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenC_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalC_addr0_0(i, j, ireuse, jreuse) );
                        }
                        findReuseFlag = true;
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr0_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                    if ( calAddrC_addr0(i, j) == calAddrC_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenC_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalC_addr0_1(i, j, ireuse, jreuse) );
                        }
                        findReuseFlag = true;
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr0_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 256; ireuse++) {
            for ( int jreuse = 0; jreuse < 256; jreuse++) {
                for ( int kreuse = 0; kreuse < 256; kreuse++) {
                        if ( calAddrC_addr0(i, j) == calAddrC_addr2(ireuse, jreuse, kreuse) ) {
                            if ( checkIntervenC_addr0_2(i, j, ireuse, jreuse, kreuse) == false) {
                              rtHistoCal(  rtCalC_addr0_2(i, j, ireuse, jreuse, kreuse) );
                            }
                            findReuseFlag = true;
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr0_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 256; ireuse++) {
            for ( int jreuse = 0; jreuse < 256; jreuse++) {
                for ( int kreuse = 0; kreuse < 256; kreuse++) {
                        if ( calAddrC_addr0(i, j) == calAddrC_addr3(ireuse, jreuse, kreuse) ) {
                            if ( checkIntervenC_addr0_3(i, j, ireuse, jreuse, kreuse) == false) {
                              rtHistoCal(  rtCalC_addr0_3(i, j, ireuse, jreuse, kreuse) );
                            }
                            findReuseFlag = true;
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                    if ( calAddrC_addr1(i, j) == calAddrC_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenC_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalC_addr1_0(i, j, ireuse, jreuse) );
                        }
                        findReuseFlag = true;
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                    if ( calAddrC_addr1(i, j) == calAddrC_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenC_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalC_addr1_1(i, j, ireuse, jreuse) );
                        }
                        findReuseFlag = true;
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 256; ireuse++) {
            for ( int jreuse = 0; jreuse < 256; jreuse++) {
                for ( int kreuse = 0; kreuse < 256; kreuse++) {
                        if ( calAddrC_addr1(i, j) == calAddrC_addr2(ireuse, jreuse, kreuse) ) {
                            if ( checkIntervenC_addr1_2(i, j, ireuse, jreuse, kreuse) == false) {
                              rtHistoCal(  rtCalC_addr1_2(i, j, ireuse, jreuse, kreuse) );
                            }
                            findReuseFlag = true;
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 256; ireuse++) {
            for ( int jreuse = 0; jreuse < 256; jreuse++) {
                for ( int kreuse = 0; kreuse < 256; kreuse++) {
                        if ( calAddrC_addr1(i, j) == calAddrC_addr3(ireuse, jreuse, kreuse) ) {
                            if ( checkIntervenC_addr1_3(i, j, ireuse, jreuse, kreuse) == false) {
                              rtHistoCal(  rtCalC_addr1_3(i, j, ireuse, jreuse, kreuse) );
                            }
                            findReuseFlag = true;
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairC_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairC_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairC_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrC_addr2(i, j, k) == calAddrC_addr2(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenC_addr2_2(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalC_addr2_2(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairC_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrC_addr2(i, j, k) == calAddrC_addr3(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenC_addr2_3(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalC_addr2_3(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairC_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairC_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairC_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrC_addr3(i, j, k) == calAddrC_addr2(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenC_addr3_2(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalC_addr3_2(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void pairC_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (256 - 0) + 0;
        int j = rand() % (256 - 0) + 0;
        int k = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (256 - 0) + 0;
            j = rand() % (256 - 0) + 0;
            k = rand() % (256 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 256; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 256; jreuse++) {
                int kreuse;
                if (ireuse == i && jreuse == j) {
                    kreuse = k + 1;
                } else {
                    kreuse = 0;
                }
                for ( ; kreuse < 256; kreuse++) {
                            if ( calAddrC_addr3(i, j, k) == calAddrC_addr3(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenC_addr3_3(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalC_addr3_3(i, j, k, ireuse, jreuse, kreuse) );
                                }
                                findReuseFlag = true;
                            }
                            if (findReuseFlag == true) {
                                break;
                            }
                        }
                        if (findReuseFlag == true) {
                            break;
                        }
                    }
                    if (findReuseFlag == true) {
                        break;
                    }
                }
            }
}
void rtDump() {
    cout << "Start to dump reuse time histogram\n";
    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {
        cout << it->first << " " << it->second << "\n";
    }
    return;
}
void dumpMR() {
    cout << "miss ratio" << endl;
    std::map<uint64_t, double>::iterator it1 = MR.begin();
    std::map<uint64_t, double>::iterator it2 = MR.begin();
    while(it1 != MR.end()) {
        while(1) {
            std::map<uint64_t, double>::iterator it3 = it2;
            ++it3;
            if (it3 == MR.end()) {
                break;
            }
            if (it1->second - it3->second < 0.00001) {
                ++it2;
            } else {
                break;
            }
        }
        cout << it1->first << " " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << " " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
int main() {
    pairA_addr0_0();
    pairA_addr0_1();
    pairA_addr1_0();
    pairA_addr1_1();
    pairC_addr0_0();
    pairC_addr0_1();
    pairC_addr0_2();
    pairC_addr0_3();
    pairC_addr1_0();
    pairC_addr1_1();
    pairC_addr1_2();
    pairC_addr1_3();
    pairC_addr2_0();
    pairC_addr2_1();
    pairC_addr2_2();
    pairC_addr2_3();
    pairC_addr3_0();
    pairC_addr3_1();
    pairC_addr3_2();
    pairC_addr3_3();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
syrk */ 
