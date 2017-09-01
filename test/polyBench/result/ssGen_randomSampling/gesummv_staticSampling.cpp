
 /* Start to analysis array index
Array index info
y.addr i
tmp.addr i
A.addr ((i * 1024) + j)
x.addr j
tmp.addr i
tmp.addr i
B.addr ((i * 1024) + j)
x.addr j
y.addr i
y.addr i
tmp.addr i
y.addr i
y.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double %alpha
double %beta
double* %A
double* %B
double* %tmp
double* %x
double* %y

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx, align 8
----array access   store double 0.000000e+00, double* %arrayidx4, align 8
----j
----Loop Bound: (0, 1024)
------array access   %9 = load double, double* %arrayidx9, align 8
------array access   %12 = load double, double* %arrayidx11, align 8
------array access   %15 = load double, double* %arrayidx14, align 8
------array access   store double %add15, double* %arrayidx17, align 8
------array access   %21 = load double, double* %arrayidx21, align 8
------array access   %24 = load double, double* %arrayidx23, align 8
------array access   %27 = load double, double* %arrayidx26, align 8
------array access   store double %add27, double* %arrayidx29, align 8
----array access   %34 = load double, double* %arrayidx31, align 8
----array access   %38 = load double, double* %arrayidx34, align 8
----array access   store double %add36, double* %arrayidx38, align 8

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
int calAddrtmp_addr0( int i) {
    int result = i;
    return result;
}
int calAddry_addr0( int i) {
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
int calAddrtmp_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddrtmp_addr2( int i, int j) {
    int result = i;
    return result;
}
int calAddrB_addr0( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrx_addr1( int i, int j) {
    int result = j;
    return result;
}
int calAddry_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddry_addr2( int i, int j) {
    int result = i;
    return result;
}
int calAddrtmp_addr3( int i) {
    int result = i;
    return result;
}
int calAddry_addr3( int i) {
    int result = i;
    return result;
}
int calAddry_addr4( int i) {
    int result = i;
    return result;
}
int rtCalA_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 0 - 0;
}
int rtCalB_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 4 - 4;
}
int rtCaltmp_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 8197 + 0 - 0;
}
int rtCaltmp_addr0_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 +  1 + 1 + (jreuse - 0) * 8 + 2;
}
int rtCaltmp_addr0_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 +  1 + 1 + (jreuse - 0) * 8 + 3;
}
int rtCaltmp_addr0_3(int i, int ireuse) {
    return (ireuse - i) * 8197 + 2 - 0;
}
int rtCaltmp_addr1_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 2;
}
int rtCaltmp_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 2 - 2;
}
int rtCaltmp_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 3 - 2;
}
int rtCaltmp_addr1_3(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 2;
}
int rtCaltmp_addr2_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 3;
}
int rtCaltmp_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 2 - 3;
}
int rtCaltmp_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 3 - 3;
}
int rtCaltmp_addr2_3(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 3;
}
int rtCaltmp_addr3_0(int i, int ireuse) {
    return (ireuse - i) * 8197 + 0 - 2;
}
int rtCaltmp_addr3_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 2;
}
int rtCaltmp_addr3_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 3;
}
int rtCaltmp_addr3_3(int i, int ireuse) {
    return (ireuse - i) * 8197 + 2 - 2;
}
int rtCalx_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 1 - 1;
}
int rtCalx_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 5 - 1;
}
int rtCalx_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 1 - 5;
}
int rtCalx_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 5 - 5;
}
int rtCaly_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 8197 + 1 - 1;
}
int rtCaly_addr0_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 +  1 + (jreuse - 0) * 8 + 6;
}
int rtCaly_addr0_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 +  1 + (jreuse - 0) * 8 + 7;
}
int rtCaly_addr0_3(int i, int ireuse) {
    return (ireuse - i) * 8197 + 3 - 1;
}
int rtCaly_addr0_4(int i, int ireuse) {
    return (ireuse - i) * 8197 + 4 - 1;
}
int rtCaly_addr1_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 6;
}
int rtCaly_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 6 - 6;
}
int rtCaly_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 7 - 6;
}
int rtCaly_addr1_3(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + 1 + (1024 - j) * 8  - 6;
}
int rtCaly_addr1_4(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + 1 + 1 + (1024 - j) * 8  - 6;
}
int rtCaly_addr2_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + (1024 - j) * 8  - 7;
}
int rtCaly_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 6 - 7;
}
int rtCaly_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - j) * 8 + 7 - 7;
}
int rtCaly_addr2_3(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + 1 + (1024 - j) * 8  - 7;
}
int rtCaly_addr2_4(int i, int j, int ireuse) {
    return (ireuse - i) * 8197 + 1 + 1 + (1024 - j) * 8  - 7;
}
int rtCaly_addr3_0(int i, int ireuse) {
    return (ireuse - i) * 8197 + 1 - 3;
}
int rtCaly_addr3_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 6;
}
int rtCaly_addr3_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 7;
}
int rtCaly_addr3_3(int i, int ireuse) {
    return (ireuse - i) * 8197 + 3 - 3;
}
int rtCaly_addr3_4(int i, int ireuse) {
    return (ireuse - i) * 8197 + 4 - 3;
}
int rtCaly_addr4_0(int i, int ireuse) {
    return (ireuse - i) * 8197 + 1 - 4;
}
int rtCaly_addr4_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 6;
}
int rtCaly_addr4_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8197 + (jreuse - 0) * 8 + 7;
}
int rtCaly_addr4_3(int i, int ireuse) {
    return (ireuse - i) * 8197 + 3 - 4;
}
int rtCaly_addr4_4(int i, int ireuse) {
    return (ireuse - i) * 8197 + 4 - 4;
}
bool checkIntervenA_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenB_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr0(i)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr0(i)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr0_3(int i, int ireuse) { 
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
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr1_3(int i, int j, int ireuse) { 
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
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
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
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr2_3(int i, int j, int ireuse) { 
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
bool checkInterventmp_addr3_0(int i, int ireuse) { 
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
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i)) {
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
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr3_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i)) {
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
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr3_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i)) {
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
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr3(iInterven) == calAddrtmp_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterventmp_addr3_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrtmp_addr0(iInterven) == calAddrtmp_addr3(i)) {
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
            if( calAddrtmp_addr1(iInterven, jInterven) == calAddrtmp_addr3(i)) {
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
            if( calAddrtmp_addr2(iInterven, jInterven) == calAddrtmp_addr3(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx_addr0_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx_addr1_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny_addr0_0(int i, int ireuse) { 
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr0_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr0_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr0_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr0_4(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr0(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr0(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr0(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr1_0(int i, int j, int ireuse) { 
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr1_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr1_3(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr1_4(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr1(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr2_0(int i, int j, int ireuse) { 
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr2_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr2_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr2_3(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr2_4(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr2(i, j)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr2(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr2(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr3_0(int i, int ireuse) { 
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr3_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr3_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr3_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr3_4(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr3(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr3(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr4_0(int i, int ireuse) { 
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr4(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr4_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr4(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr4_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr4(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr4_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr4(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddry_addr4(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    return false;
}
bool checkInterveny_addr4_4(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr0(iInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr1(iInterven, jInterven) == calAddry_addr4(i)) {
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
            if( calAddry_addr2(iInterven, jInterven) == calAddry_addr4(i)) {
                return true;
            }
        }
    }
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddry_addr3(iInterven) == calAddry_addr4(i)) {
            return true;
        }
    }
    return false;
}
void pairA_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairB_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrB_addr0(i, j) == calAddrB_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenB_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalB_addr0_0(i, j, ireuse, jreuse) );
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
void pairtmp_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
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
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
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
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
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
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            if ( calAddrtmp_addr0(i) == calAddrtmp_addr3(ireuse) ) {
                if ( checkInterventmp_addr0_3(i, ireuse) == false) {
                  rtHistoCal(  rtCaltmp_addr0_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairtmp_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddrtmp_addr1(i, j) == calAddrtmp_addr3(ireuse) ) {
                    if ( checkInterventmp_addr1_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaltmp_addr1_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairtmp_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairtmp_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddrtmp_addr2(i, j) == calAddrtmp_addr3(ireuse) ) {
                    if ( checkInterventmp_addr2_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaltmp_addr2_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairtmp_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrtmp_addr3(i) == calAddrtmp_addr0(ireuse) ) {
                if ( checkInterventmp_addr3_0(i, ireuse) == false) {
                  rtHistoCal(  rtCaltmp_addr3_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairtmp_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrtmp_addr3(i) == calAddrtmp_addr1(ireuse, jreuse) ) {
                    if ( checkInterventmp_addr3_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaltmp_addr3_1(i, ireuse, jreuse) );
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
void pairtmp_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrtmp_addr3(i) == calAddrtmp_addr2(ireuse, jreuse) ) {
                    if ( checkInterventmp_addr3_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaltmp_addr3_2(i, ireuse, jreuse) );
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
void pairtmp_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrtmp_addr3(i) == calAddrtmp_addr3(ireuse) ) {
                if ( checkInterventmp_addr3_3(i, ireuse) == false) {
                  rtHistoCal(  rtCaltmp_addr3_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairx_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairx_addr0_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr0(i, j) == calAddrx_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr0_1(i, j, ireuse, jreuse) );
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
void pairx_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr1(i, j) == calAddrx_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr1_0(i, j, ireuse, jreuse) );
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
void pairx_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 0;
            }
            for ( ; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr1(i, j) == calAddrx_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr1_1(i, j, ireuse, jreuse) );
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
void pairy_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
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
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
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
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
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
void pairy_addr0_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            if ( calAddry_addr0(i) == calAddry_addr3(ireuse) ) {
                if ( checkInterveny_addr0_3(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr0_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr0_4() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            if ( calAddry_addr0(i) == calAddry_addr4(ireuse) ) {
                if ( checkInterveny_addr0_4(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr0_4(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                if ( calAddry_addr1(i, j) == calAddry_addr0(ireuse) ) {
                    if ( checkInterveny_addr1_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr1_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairy_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairy_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddry_addr1(i, j) == calAddry_addr3(ireuse) ) {
                    if ( checkInterveny_addr1_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr1_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr1_4() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddry_addr1(i, j) == calAddry_addr4(ireuse) ) {
                    if ( checkInterveny_addr1_4(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr1_4(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
                if ( calAddry_addr2(i, j) == calAddry_addr0(ireuse) ) {
                    if ( checkInterveny_addr2_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr2_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairy_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
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
void pairy_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddry_addr2(i, j) == calAddry_addr3(ireuse) ) {
                    if ( checkInterveny_addr2_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr2_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr2_4() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        int j = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            j = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
                if ( calAddry_addr2(i, j) == calAddry_addr4(ireuse) ) {
                    if ( checkInterveny_addr2_4(i, j, ireuse) == false) {
                      rtHistoCal(  rtCaly_addr2_4(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairy_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr3(i) == calAddry_addr0(ireuse) ) {
                if ( checkInterveny_addr3_0(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr3_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr3(i) == calAddry_addr1(ireuse, jreuse) ) {
                    if ( checkInterveny_addr3_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr3_1(i, ireuse, jreuse) );
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
void pairy_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr3(i) == calAddry_addr2(ireuse, jreuse) ) {
                    if ( checkInterveny_addr3_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr3_2(i, ireuse, jreuse) );
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
void pairy_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr3(i) == calAddry_addr3(ireuse) ) {
                if ( checkInterveny_addr3_3(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr3_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr3_4() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            if ( calAddry_addr3(i) == calAddry_addr4(ireuse) ) {
                if ( checkInterveny_addr3_4(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr3_4(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr4_0() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr4(i) == calAddry_addr0(ireuse) ) {
                if ( checkInterveny_addr4_0(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr4_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr4_1() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr4(i) == calAddry_addr1(ireuse, jreuse) ) {
                    if ( checkInterveny_addr4_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr4_1(i, ireuse, jreuse) );
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
void pairy_addr4_2() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddry_addr4(i) == calAddry_addr2(ireuse, jreuse) ) {
                    if ( checkInterveny_addr4_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCaly_addr4_2(i, ireuse, jreuse) );
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
void pairy_addr4_3() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr4(i) == calAddry_addr3(ireuse) ) {
                if ( checkInterveny_addr4_3(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr4_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairy_addr4_4() {
    set<string> record;
    for ( int s = 0; s < 10000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddry_addr4(i) == calAddry_addr4(ireuse) ) {
                if ( checkInterveny_addr4_4(i, ireuse) == false) {
                  rtHistoCal(  rtCaly_addr4_4(i, ireuse) );
                }
                findReuseFlag = true;
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
    pairB_addr0_0();
    pairtmp_addr0_0();
    pairtmp_addr0_1();
    pairtmp_addr0_2();
    pairtmp_addr0_3();
    pairtmp_addr1_0();
    pairtmp_addr1_1();
    pairtmp_addr1_2();
    pairtmp_addr1_3();
    pairtmp_addr2_0();
    pairtmp_addr2_1();
    pairtmp_addr2_2();
    pairtmp_addr2_3();
    pairtmp_addr3_0();
    pairtmp_addr3_1();
    pairtmp_addr3_2();
    pairtmp_addr3_3();
    pairx_addr0_0();
    pairx_addr0_1();
    pairx_addr1_0();
    pairx_addr1_1();
    pairy_addr0_0();
    pairy_addr0_1();
    pairy_addr0_2();
    pairy_addr0_3();
    pairy_addr0_4();
    pairy_addr1_0();
    pairy_addr1_1();
    pairy_addr1_2();
    pairy_addr1_3();
    pairy_addr1_4();
    pairy_addr2_0();
    pairy_addr2_1();
    pairy_addr2_2();
    pairy_addr2_3();
    pairy_addr2_4();
    pairy_addr3_0();
    pairy_addr3_1();
    pairy_addr3_2();
    pairy_addr3_3();
    pairy_addr3_4();
    pairy_addr4_0();
    pairy_addr4_1();
    pairy_addr4_2();
    pairy_addr4_3();
    pairy_addr4_4();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
gesummv */ 
