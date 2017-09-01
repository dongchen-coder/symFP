
 /* Start to analysis array index
Array index info
B.addr ((k * 256) + j)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)
A.addr ((i * 256) + k)
C.addr ((i * 256) + j)
C.addr ((i * 256) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
double %alpha
double %beta
double* %A
double* %B
double* %C

 Finish to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 256)
----j
----Loop Bound: (0, 256)
------array access   %6 = load double, double* %arrayidx, align 8
------array access   store double %mul4, double* %arrayidx, align 8
------k
------Loop Bound: (0, 256)
--------array access   %12 = load double, double* %arrayidx11, align 8
--------array access   %16 = load double, double* %arrayidx16, align 8
--------array access   %20 = load double, double* %arrayidx21, align 8
--------array access   store double %add22, double* %arrayidx21, align 8

Finish analysis loops */ 
 // Start to generating Static Sampling Code
#include <map>
#include <iostream>
using namespace std;
std::map<uint64_t, uint64_t> rtHisto;
void rtHistoCal( int rt) {
    if (rtHisto.find(rt) == rtHisto.end()) { 
        rtHisto[rt] = 1;
    } else {
        rtHisto[rt] += 1;
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
int calAddrB_addr0( int i, int j, int k) {
    int result = ((k * 256) + j);
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
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 0 - 0;
}
int rtCalB_addr0_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 1 - 1;
}
int rtCalC_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + 0 - 0;
}
int rtCalC_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + 1 - 0;
}
int rtCalC_addr0_2(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 +  1 + 1 + (kreuse - 0) * 4 + 2;
}
int rtCalC_addr0_3(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 +  1 + 1 + (kreuse - 0) * 4 + 3;
}
int rtCalC_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + 0 - 1;
}
int rtCalC_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + 1 - 1;
}
int rtCalC_addr1_2(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 +  1 + (kreuse - 0) * 4 + 2;
}
int rtCalC_addr1_3(int i, int j, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 +  1 + (kreuse - 0) * 4 + 3;
}
int rtCalC_addr2_0(int i, int j, int k, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (256 - k) * 4  - 2;
}
int rtCalC_addr2_1(int i, int j, int k, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (256 - k) * 4  - 2;
}
int rtCalC_addr2_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 2 - 2;
}
int rtCalC_addr2_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 3 - 2;
}
int rtCalC_addr3_0(int i, int j, int k, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (256 - k) * 4  - 3;
}
int rtCalC_addr3_1(int i, int j, int k, int ireuse, int jreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (256 - k) * 4  - 3;
}
int rtCalC_addr3_2(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 2 - 3;
}
int rtCalC_addr3_3(int i, int j, int k, int ireuse, int jreuse, int kreuse) {
    return (ireuse - i) * 262656 + (jreuse - j) * 1026 + (kreuse - k) * 4 + 3 - 3;
}
bool checkIntervenA_addr0_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
    return false;
}
bool checkIntervenB_addr0_0(int i, int j, int k, int ireuse, int jreuse, int kreuse) { 
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
bool checkIntervenC_addr0_2(int i, int j, int ireuse, int jreuse, int kreuse) { 
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr0(i, j)) {
                return true;
            }
        }
    }
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = 0;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
bool checkIntervenC_addr1_2(int i, int j, int ireuse, int jreuse, int kreuse) { 
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr1(i, j)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
bool checkIntervenC_addr2_1(int i, int j, int k, int ireuse, int jreuse) { 
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr2(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
bool checkIntervenC_addr3_1(int i, int j, int k, int ireuse, int jreuse) { 
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            int kInterven;
            int kIntervenUB;
            if (iInterven == i && jInterven == j) {
                kInterven = k + 1 ;
            } else {
                kInterven = 0;
            }
            if (iInterven == ireuse && jInterven == jreuse - 1 ) {
                kIntervenUB = 256;
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr0(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
            if( calAddrC_addr1(iInterven, jInterven) == calAddrC_addr3(i, j, k)) {
                return true;
            }
        }
    }
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
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
    }
}
void pairB_addr0_0() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
                            if ( calAddrB_addr0(i, j, k) == calAddrB_addr0(ireuse, jreuse, kreuse) ) {
                                if ( checkIntervenB_addr0_0(i, j, k, ireuse, jreuse, kreuse) == false) {
                                  rtHistoCal(  rtCalB_addr0_0(i, j, k, ireuse, jreuse, kreuse) );
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
    }
}
void pairC_addr0_0() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
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
}
void pairC_addr0_1() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
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
}
void pairC_addr0_2() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 256; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 256; jreuse++) {
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
}
void pairC_addr0_3() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 256; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 256; jreuse++) {
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
}
void pairC_addr1_0() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
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
}
void pairC_addr1_1() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
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
}
void pairC_addr1_2() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 256; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 256; jreuse++) {
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
}
void pairC_addr1_3() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 256; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 256; jreuse++) {
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
}
void pairC_addr2_0() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
                bool findReuseFlag = false;
                for ( int ireuse = i; ireuse < 256; ireuse++) {
                    int jreuse;
                    if (ireuse == i) {
                        jreuse = j + 1;
                    } else {
                        jreuse = 0;
                    }
                    for ( ; jreuse < 256; jreuse++) {
                        if ( calAddrC_addr2(i, j, k) == calAddrC_addr0(ireuse, jreuse) ) {
                            if ( checkIntervenC_addr2_0(i, j, k, ireuse, jreuse) == false) {
                              rtHistoCal(  rtCalC_addr2_0(i, j, k, ireuse, jreuse) );
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
    }
}
void pairC_addr2_1() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
                bool findReuseFlag = false;
                for ( int ireuse = i; ireuse < 256; ireuse++) {
                    int jreuse;
                    if (ireuse == i) {
                        jreuse = j + 1;
                    } else {
                        jreuse = 0;
                    }
                    for ( ; jreuse < 256; jreuse++) {
                        if ( calAddrC_addr2(i, j, k) == calAddrC_addr1(ireuse, jreuse) ) {
                            if ( checkIntervenC_addr2_1(i, j, k, ireuse, jreuse) == false) {
                              rtHistoCal(  rtCalC_addr2_1(i, j, k, ireuse, jreuse) );
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
    }
}
void pairC_addr2_2() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
    }
}
void pairC_addr2_3() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
    }
}
void pairC_addr3_0() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
                bool findReuseFlag = false;
                for ( int ireuse = i; ireuse < 256; ireuse++) {
                    int jreuse;
                    if (ireuse == i) {
                        jreuse = j + 1;
                    } else {
                        jreuse = 0;
                    }
                    for ( ; jreuse < 256; jreuse++) {
                        if ( calAddrC_addr3(i, j, k) == calAddrC_addr0(ireuse, jreuse) ) {
                            if ( checkIntervenC_addr3_0(i, j, k, ireuse, jreuse) == false) {
                              rtHistoCal(  rtCalC_addr3_0(i, j, k, ireuse, jreuse) );
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
    }
}
void pairC_addr3_1() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
                bool findReuseFlag = false;
                for ( int ireuse = i; ireuse < 256; ireuse++) {
                    int jreuse;
                    if (ireuse == i) {
                        jreuse = j + 1;
                    } else {
                        jreuse = 0;
                    }
                    for ( ; jreuse < 256; jreuse++) {
                        if ( calAddrC_addr3(i, j, k) == calAddrC_addr1(ireuse, jreuse) ) {
                            if ( checkIntervenC_addr3_1(i, j, k, ireuse, jreuse) == false) {
                              rtHistoCal(  rtCalC_addr3_1(i, j, k, ireuse, jreuse) );
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
    }
}
void pairC_addr3_2() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
    }
}
void pairC_addr3_3() {
    for ( int i = 0; i < 256; i += 1 /0.010000) {
        for ( int j = 0; j < 256; j += 1 /0.010000) {
            for ( int k = 0; k < 256; k += 1 /0.010000) {
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
    }
}
void rtDump() {
    cout << "Start to dump reuse time histogram\n";
    for (map<uint64_t, uint64_t>::iterator it = rtHisto.begin(), eit = rtHisto.end(); it != eit; ++it) {
        cout << it->first << " " << it->second << "\n";
    }
    return;
}
int main() {
    cout << " check pair A_addr0 A_addr0\n ";
    pairA_addr0_0();
    cout << " check pair B_addr0 B_addr0\n ";
    pairB_addr0_0();
    cout << " check pair C_addr0 C_addr0\n ";
    pairC_addr0_0();
    cout << " check pair C_addr0 C_addr1\n ";
    pairC_addr0_1();
    cout << " check pair C_addr0 C_addr2\n ";
    pairC_addr0_2();
    cout << " check pair C_addr0 C_addr3\n ";
    pairC_addr0_3();
    cout << " check pair C_addr1 C_addr0\n ";
    pairC_addr1_0();
    cout << " check pair C_addr1 C_addr1\n ";
    pairC_addr1_1();
    cout << " check pair C_addr1 C_addr2\n ";
    pairC_addr1_2();
    cout << " check pair C_addr1 C_addr3\n ";
    pairC_addr1_3();
    cout << " check pair C_addr2 C_addr0\n ";
    pairC_addr2_0();
    cout << " check pair C_addr2 C_addr1\n ";
    pairC_addr2_1();
    cout << " check pair C_addr2 C_addr2\n ";
    pairC_addr2_2();
    cout << " check pair C_addr2 C_addr3\n ";
    pairC_addr2_3();
    cout << " check pair C_addr3 C_addr0\n ";
    pairC_addr3_0();
    cout << " check pair C_addr3 C_addr1\n ";
    pairC_addr3_1();
    cout << " check pair C_addr3 C_addr2\n ";
    pairC_addr3_2();
    cout << " check pair C_addr3 C_addr3\n ";
    pairC_addr3_3();
    rtDump();
    return 0;
}
 /* Start to analyze function:  
gemm */ 
