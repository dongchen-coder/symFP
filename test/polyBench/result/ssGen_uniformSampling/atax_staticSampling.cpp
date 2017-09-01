
 /* Start to analysis array index
Array index info
tmp.addr i
y.addr i
tmp.addr i
A.addr ((i * 1024) + j)
x.addr j
tmp.addr i
y.addr j
A.addr ((i * 1024) + j)
tmp.addr i
y.addr j

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %x
double* %y
double* %tmp

 Finish to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx, align 8
--i
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx7, align 8
----j
----Loop Bound: (0, 1024)
------array access   %10 = load double, double* %arrayidx12, align 8
------array access   %14 = load double, double* %arrayidx14, align 8
------array access   %17 = load double, double* %arrayidx16, align 8
------array access   store double %add18, double* %arrayidx20, align 8
----j
----Loop Bound: (0, 1024)
------array access   %24 = load double, double* %arrayidx28, align 8
------array access   %28 = load double, double* %arrayidx32, align 8
------array access   %31 = load double, double* %arrayidx34, align 8
------array access   store double %add36, double* %arrayidx38, align 8

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
int calAddry_addr0( int i) {
    int result = i;
    return result;
}
int calAddrtmp_addr0( int i) {
    int result = i;
    return result;
}
int calAddrtmp_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddrA_addr0( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrx_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrtmp_addr2( int i, int j) {
    int result = i;
    return result;
}
int calAddry_addr1( int i, int j) {
    int result = j;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrtmp_addr3( int i, int j) {
    int result = i;
    return result;
}
int calAddry_addr2( int i, int j) {
    int result = j;
    return result;
}
int rtCalA_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 1 - 1;
}
int rtCalA_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 1 - 1;
}
int rtCalA_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 1 - 1;
}
int rtCalA_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 1 - 1;
}
int rtCaltmp_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 8193 + 0 - 0;
}
int rtCaltmp_addr0_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 +  1 + (jreuse - 0) * 4 + 0;
}
int rtCaltmp_addr0_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 +  1 + (jreuse - 0) * 4 + 3;
}
int rtCaltmp_addr0_3(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 +  1 + 4 + (jreuse - 0) * 4 + 2;
}
int rtCaltmp_addr1_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4  - 0;
}
int rtCaltmp_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 0 - 0;
}
int rtCaltmp_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 3 - 0;
}
int rtCaltmp_addr1_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 2 - 0;
}
int rtCaltmp_addr2_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4  - 3;
}
int rtCaltmp_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 0 - 3;
}
int rtCaltmp_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 3 - 3;
}
int rtCaltmp_addr2_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 2 - 3;
}
int rtCaltmp_addr3_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4  - 2;
}
int rtCaltmp_addr3_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 0 - 2;
}
int rtCaltmp_addr3_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 4 + (jreuse - 0) * 4 + 3 - 2;
}
int rtCaltmp_addr3_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 2 - 2;
}
int rtCalx_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 2 - 2;
}
int rtCaly_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 1 + 0 - 0;
}
int rtCaly_addr0_1(int i, int ireuse, int jreuse) {
    return (1024 - i) * 1 + (ireuse - 0) * 8193 + (jreuse - 0) * 4 + 0 - 0;
}
int rtCaly_addr0_2(int i, int ireuse, int jreuse) {
    return (1024 - i) * 1 + (ireuse - 0) * 8193 + (jreuse - 0) * 4 + 3 - 0;
}
int rtCaly_addr1_0(int i, int j, int ireuse) {
    return (1024 - i) * 8193 + (ireuse - 0) * 1 + (1024 - j) * 4 + 0 - 0;
}
int rtCaly_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 0 - 0;
}
int rtCaly_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 3 - 0;
}
int rtCaly_addr2_0(int i, int j, int ireuse) {
    return (1024 - i) * 8193 + (ireuse - 0) * 1 + (1024 - j) * 4 + 0 - 3;
}
int rtCaly_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 0 - 3;
}
int rtCaly_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 4 + 3 - 3;
}
bool checkIntervenA_addr0_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr1(iInterven, jInterven) == calAddrA_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr0_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr0(iInterven, jInterven) == calAddrA_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr1_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr1(iInterven, jInterven) == calAddrA_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr0(iInterven, jInterven) == calAddrA_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr0_0(int i, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr0_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr0_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr0_3(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr1_0(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr1_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr1_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
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
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
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
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr2_0(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr2_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
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
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr2_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr2_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
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
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
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
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr3_0(int i, int j, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr3_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i, j)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr3_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i, j)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse - 1 ) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr3(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterventmp_addr3_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i, j)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i + 1 ) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == ireuse) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkInterveny_addr0_0(int i, int ireuse) { 
    return false;
}
bool checkInterveny_addr0_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
            return true;
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
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr0_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
            return true;
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr1_0(int i, int j, int ireuse) { 
    return false;
}
bool checkInterveny_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr1_2(int i, int j, int ireuse, int jreuse) { 
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr2_0(int i, int j, int ireuse) { 
    return false;
}
bool checkInterveny_addr2_1(int i, int j, int ireuse, int jreuse) { 
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr2_2(int i, int j, int ireuse, int jreuse) { 
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
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
void pairA_addr0_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr0(i, j) == calAddrA_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr0_0(i, j, ireuse, jreuse) );
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
void pairA_addr0_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr0(i, j) == calAddrA_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr0_1(i, j, ireuse, jreuse) );
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
void pairA_addr1_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr1(i, j) == calAddrA_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr1_0(i, j, ireuse, jreuse) );
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
void pairA_addr1_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr1(i, j) == calAddrA_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr1_1(i, j, ireuse, jreuse) );
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
void pairtmp_addr0_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrtmp_addr0(i) == calAddrtmp_addr0(ireuse) ) {
                if ( checkInterventmp_addr0_0(i, ireuse) == false) {
                  rtHistoCal(  rtCaltmp_addr0_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairtmp_addr0_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrtmp_addr0(i) == calAddrtmp_addr1(ireuse, jreuse) ) {
                    if ( checkInterventmp_addr0_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaltmp_addr0_1(i, ireuse, jreuse) );
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
void pairtmp_addr0_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrtmp_addr0(i) == calAddrtmp_addr2(ireuse, jreuse) ) {
                    if ( checkInterventmp_addr0_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaltmp_addr0_2(i, ireuse, jreuse) );
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
void pairtmp_addr0_3() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrtmp_addr0(i) == calAddrtmp_addr3(ireuse, jreuse) ) {
                    if ( checkInterventmp_addr0_3(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaltmp_addr0_3(i, ireuse, jreuse) );
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
void pairtmp_addr1_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                if ( calAddrtmp_addr1(i, j) == calAddrtmp_addr0(ireuse) ) {
                    if ( checkInterventmp_addr1_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaltmp_addr1_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
    }
}
void pairtmp_addr1_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr1(i, j) == calAddrtmp_addr1(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr1_1(i, j, ireuse, jreuse) );
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
void pairtmp_addr1_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr1(i, j) == calAddrtmp_addr2(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr1_2(i, j, ireuse, jreuse) );
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
void pairtmp_addr1_3() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr1(i, j) == calAddrtmp_addr3(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr1_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr1_3(i, j, ireuse, jreuse) );
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
void pairtmp_addr2_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                if ( calAddrtmp_addr2(i, j) == calAddrtmp_addr0(ireuse) ) {
                    if ( checkInterventmp_addr2_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaltmp_addr2_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
    }
}
void pairtmp_addr2_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr2(i, j) == calAddrtmp_addr1(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr2_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr2_1(i, j, ireuse, jreuse) );
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
void pairtmp_addr2_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr2(i, j) == calAddrtmp_addr2(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr2_2(i, j, ireuse, jreuse) );
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
void pairtmp_addr2_3() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr2(i, j) == calAddrtmp_addr3(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr2_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr2_3(i, j, ireuse, jreuse) );
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
void pairtmp_addr3_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                if ( calAddrtmp_addr3(i, j) == calAddrtmp_addr0(ireuse) ) {
                    if ( checkInterventmp_addr3_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaltmp_addr3_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
    }
}
void pairtmp_addr3_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr3(i, j) == calAddrtmp_addr1(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr3_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr3_1(i, j, ireuse, jreuse) );
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
void pairtmp_addr3_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr3(i, j) == calAddrtmp_addr2(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr3_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr3_2(i, j, ireuse, jreuse) );
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
void pairtmp_addr3_3() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrtmp_addr3(i, j) == calAddrtmp_addr3(ireuse, jreuse) ) {
                        if ( checkInterventmp_addr3_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaltmp_addr3_3(i, j, ireuse, jreuse) );
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
void pairx_addr0_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr0(i, j) == calAddrx_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr0_0(i, j, ireuse, jreuse) );
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
void pairy_addr0_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr0(i) == calAddry_addr0(ireuse) ) {
                if ( checkInterveny_addr0_0(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr0_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr0_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr0(i) == calAddry_addr1(ireuse, jreuse) ) {
                    if ( checkInterveny_addr0_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr0_1(i, ireuse, jreuse) );
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
void pairy_addr0_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        bool findReuseFlag = false;
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr0(i) == calAddry_addr2(ireuse, jreuse) ) {
                    if ( checkInterveny_addr0_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr0_2(i, ireuse, jreuse) );
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
void pairy_addr1_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
    break;
}
    break;
}
}
void pairy_addr1_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddry_addr1(i, j) == calAddry_addr1(ireuse, jreuse) ) {
                        if ( checkInterveny_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly_addr1_1(i, j, ireuse, jreuse) );
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
void pairy_addr1_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddry_addr1(i, j) == calAddry_addr2(ireuse, jreuse) ) {
                        if ( checkInterveny_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly_addr1_2(i, j, ireuse, jreuse) );
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
void pairy_addr2_0() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
    break;
}
    break;
}
}
void pairy_addr2_1() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddry_addr2(i, j) == calAddry_addr1(ireuse, jreuse) ) {
                        if ( checkInterveny_addr2_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly_addr2_1(i, j, ireuse, jreuse) );
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
void pairy_addr2_2() {
    for ( int i = 0; i < 1024; i += 1 /0.010000) {
        for ( int j = 0; j < 1024; j += 1 /0.010000) {
            bool findReuseFlag = false;
            for ( int ireuse = i; ireuse < 1024; ireuse++) {
                int jreuse;
                if (ireuse == i) {
                    jreuse = j + 1;
                } else {
                    jreuse = 0;
                }
                for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddry_addr2(i, j) == calAddry_addr2(ireuse, jreuse) ) {
                        if ( checkInterveny_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly_addr2_2(i, j, ireuse, jreuse) );
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
    cout << " check pair A_addr0 A_addr1\n ";
    pairA_addr0_1();
    cout << " check pair A_addr1 A_addr0\n ";
    pairA_addr1_0();
    cout << " check pair A_addr1 A_addr1\n ";
    pairA_addr1_1();
    cout << " check pair tmp_addr0 tmp_addr0\n ";
    pairtmp_addr0_0();
    cout << " check pair tmp_addr0 tmp_addr1\n ";
    pairtmp_addr0_1();
    cout << " check pair tmp_addr0 tmp_addr2\n ";
    pairtmp_addr0_2();
    cout << " check pair tmp_addr0 tmp_addr3\n ";
    pairtmp_addr0_3();
    cout << " check pair tmp_addr1 tmp_addr0\n ";
    pairtmp_addr1_0();
    cout << " check pair tmp_addr1 tmp_addr1\n ";
    pairtmp_addr1_1();
    cout << " check pair tmp_addr1 tmp_addr2\n ";
    pairtmp_addr1_2();
    cout << " check pair tmp_addr1 tmp_addr3\n ";
    pairtmp_addr1_3();
    cout << " check pair tmp_addr2 tmp_addr0\n ";
    pairtmp_addr2_0();
    cout << " check pair tmp_addr2 tmp_addr1\n ";
    pairtmp_addr2_1();
    cout << " check pair tmp_addr2 tmp_addr2\n ";
    pairtmp_addr2_2();
    cout << " check pair tmp_addr2 tmp_addr3\n ";
    pairtmp_addr2_3();
    cout << " check pair tmp_addr3 tmp_addr0\n ";
    pairtmp_addr3_0();
    cout << " check pair tmp_addr3 tmp_addr1\n ";
    pairtmp_addr3_1();
    cout << " check pair tmp_addr3 tmp_addr2\n ";
    pairtmp_addr3_2();
    cout << " check pair tmp_addr3 tmp_addr3\n ";
    pairtmp_addr3_3();
    cout << " check pair x_addr0 x_addr0\n ";
    pairx_addr0_0();
    cout << " check pair y_addr0 y_addr0\n ";
    pairy_addr0_0();
    cout << " check pair y_addr0 y_addr1\n ";
    pairy_addr0_1();
    cout << " check pair y_addr0 y_addr2\n ";
    pairy_addr0_2();
    cout << " check pair y_addr1 y_addr0\n ";
    pairy_addr1_0();
    cout << " check pair y_addr1 y_addr1\n ";
    pairy_addr1_1();
    cout << " check pair y_addr1 y_addr2\n ";
    pairy_addr1_2();
    cout << " check pair y_addr2 y_addr0\n ";
    pairy_addr2_0();
    cout << " check pair y_addr2 y_addr1\n ";
    pairy_addr2_1();
    cout << " check pair y_addr2 y_addr2\n ";
    pairy_addr2_2();
    rtDump();
    return 0;
}
 /* Start to analyze function:  
atax_cpu */ 
