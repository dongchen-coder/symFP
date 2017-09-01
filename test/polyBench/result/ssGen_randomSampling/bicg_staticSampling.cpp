
 /* Start to analysis array index
Array index info
q.addr i
s.addr i
s.addr j
r.addr i
A.addr ((i * 1024) + j)
s.addr j
q.addr i
A.addr ((i * 1024) + j)
p.addr j
q.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %r
double* %s
double* %p
double* %q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx, align 8
--i
--Loop Bound: (0, 1024)
----array access   store double 0.000000e+00, double* %arrayidx5, align 8
----j
----Loop Bound: (0, 1024)
------array access   %10 = load double, double* %arrayidx10, align 8
------array access   %13 = load double, double* %arrayidx12, align 8
------array access   %17 = load double, double* %arrayidx14, align 8
------array access   store double %add16, double* %arrayidx18, align 8
------array access   %22 = load double, double* %arrayidx20, align 8
------array access   %26 = load double, double* %arrayidx24, align 8
------array access   %29 = load double, double* %arrayidx26, align 8
------array access   store double %add28, double* %arrayidx30, align 8

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
int calAddrs_addr0( int i) {
    int result = i;
    return result;
}
int calAddrq_addr0( int i) {
    int result = i;
    return result;
}
int calAddrs_addr1( int i, int j) {
    int result = j;
    return result;
}
int calAddrr_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddrA_addr0( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrs_addr2( int i, int j) {
    int result = j;
    return result;
}
int calAddrq_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrp_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrq_addr2( int i, int j) {
    int result = i;
    return result;
}
int rtCalA_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 2 - 2;
}
int rtCalA_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 5 - 2;
}
int rtCalA_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 2 - 5;
}
int rtCalA_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 5 - 5;
}
int rtCalp_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 6 - 6;
}
int rtCalq_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 8193 + 0 - 0;
}
int rtCalq_addr0_1(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 +  1 + (jreuse - 0) * 8 + 4;
}
int rtCalq_addr0_2(int i, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 +  1 + (jreuse - 0) * 8 + 7;
}
int rtCalq_addr1_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 8  - 4;
}
int rtCalq_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 4 - 4;
}
int rtCalq_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 7 - 4;
}
int rtCalq_addr2_0(int i, int j, int ireuse) {
    return (ireuse - i) * 8193 + (1024 - j) * 8  - 7;
}
int rtCalq_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 4 - 7;
}
int rtCalq_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 7 - 7;
}
int rtCalr_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 1 - 1;
}
int rtCals_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 1 + 0 - 0;
}
int rtCals_addr0_1(int i, int ireuse, int jreuse) {
    return (1024 - i) * 1 + (ireuse - 0) * 8193 + (jreuse - 0) * 8 + 0 - 0;
}
int rtCals_addr0_2(int i, int ireuse, int jreuse) {
    return (1024 - i) * 1 + (ireuse - 0) * 8193 + (jreuse - 0) * 8 + 3 - 0;
}
int rtCals_addr1_0(int i, int j, int ireuse) {
    return (1024 - i) * 8193 + (ireuse - 0) * 1 + (1024 - j) * 8 + 0 - 0;
}
int rtCals_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 0 - 0;
}
int rtCals_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 3 - 0;
}
int rtCals_addr2_0(int i, int j, int ireuse) {
    return (1024 - i) * 8193 + (ireuse - 0) * 1 + (1024 - j) * 8 + 0 - 3;
}
int rtCals_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 0 - 3;
}
int rtCals_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 8193 + (jreuse - j) * 8 + 3 - 3;
}
bool checkIntervenA_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
            jIntervenUB = jreuse;
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
            if( calAddrA_addr1(iInterven, jInterven) == calAddrA_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrA_addr0(iInterven, jInterven) == calAddrA_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenp_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenq_addr0_0(int i, int ireuse) { 
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr0(i)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr0_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr0(i)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr0_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr0(i)) {
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr1_0(int i, int j, int ireuse) { 
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr1(i, j)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr1(i, j)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr1_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr1(i, j)) {
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr2_0(int i, int j, int ireuse) { 
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr2(i, j)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr2_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr2(i, j)) {
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
            if( calAddrq_addr2(iInterven, jInterven) == calAddrq_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenq_addr2_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        if( calAddrq_addr0(iInterven) == calAddrq_addr2(i, j)) {
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
            if( calAddrq_addr1(iInterven, jInterven) == calAddrq_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenr_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervens_addr0_0(int i, int ireuse) { 
    return false;
}
bool checkIntervens_addr0_1(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddrs_addr0(iInterven) == calAddrs_addr0(i)) {
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
            if( calAddrs_addr2(iInterven, jInterven) == calAddrs_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervens_addr0_2(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddrs_addr0(iInterven) == calAddrs_addr0(i)) {
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
            if( calAddrs_addr1(iInterven, jInterven) == calAddrs_addr0(i)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervens_addr1_0(int i, int j, int ireuse) { 
    return false;
}
bool checkIntervens_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrs_addr2(iInterven, jInterven) == calAddrs_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervens_addr1_2(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrs_addr1(iInterven, jInterven) == calAddrs_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervens_addr2_0(int i, int j, int ireuse) { 
    return false;
}
bool checkIntervens_addr2_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrs_addr2(iInterven, jInterven) == calAddrs_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervens_addr2_2(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrs_addr1(iInterven, jInterven) == calAddrs_addr2(i, j)) {
                return true;
            }
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
void pairA_addr0_1() {
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
void pairA_addr1_0() {
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
void pairA_addr1_1() {
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
void pairp_addr0_0() {
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
                    if ( calAddrp_addr0(i, j) == calAddrp_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenp_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalp_addr0_0(i, j, ireuse, jreuse) );
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
void pairq_addr0_0() {
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
            if ( calAddrq_addr0(i) == calAddrq_addr0(ireuse) ) {
                if ( checkIntervenq_addr0_0(i, ireuse) == false) {
                  rtHistoCal(  rtCalq_addr0_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairq_addr0_1() {
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
                if ( calAddrq_addr0(i) == calAddrq_addr1(ireuse, jreuse) ) {
                    if ( checkIntervenq_addr0_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCalq_addr0_1(i, ireuse, jreuse) );
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
void pairq_addr0_2() {
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
                if ( calAddrq_addr0(i) == calAddrq_addr2(ireuse, jreuse) ) {
                    if ( checkIntervenq_addr0_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCalq_addr0_2(i, ireuse, jreuse) );
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
void pairq_addr1_0() {
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
                if ( calAddrq_addr1(i, j) == calAddrq_addr0(ireuse) ) {
                    if ( checkIntervenq_addr1_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalq_addr1_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairq_addr1_1() {
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
                    if ( calAddrq_addr1(i, j) == calAddrq_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenq_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalq_addr1_1(i, j, ireuse, jreuse) );
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
void pairq_addr1_2() {
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
                    if ( calAddrq_addr1(i, j) == calAddrq_addr2(ireuse, jreuse) ) {
                        if ( checkIntervenq_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalq_addr1_2(i, j, ireuse, jreuse) );
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
void pairq_addr2_0() {
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
                if ( calAddrq_addr2(i, j) == calAddrq_addr0(ireuse) ) {
                    if ( checkIntervenq_addr2_0(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalq_addr2_0(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairq_addr2_1() {
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
                    if ( calAddrq_addr2(i, j) == calAddrq_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenq_addr2_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalq_addr2_1(i, j, ireuse, jreuse) );
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
void pairq_addr2_2() {
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
                    if ( calAddrq_addr2(i, j) == calAddrq_addr2(ireuse, jreuse) ) {
                        if ( checkIntervenq_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalq_addr2_2(i, j, ireuse, jreuse) );
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
void pairr_addr0_0() {
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
                    if ( calAddrr_addr0(i, j) == calAddrr_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenr_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalr_addr0_0(i, j, ireuse, jreuse) );
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
void pairs_addr0_0() {
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
            if ( calAddrs_addr0(i) == calAddrs_addr0(ireuse) ) {
                if ( checkIntervens_addr0_0(i, ireuse) == false) {
                  rtHistoCal(  rtCals_addr0_0(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairs_addr0_1() {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrs_addr0(i) == calAddrs_addr1(ireuse, jreuse) ) {
                    if ( checkIntervens_addr0_1(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCals_addr0_1(i, ireuse, jreuse) );
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
void pairs_addr0_2() {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                if ( calAddrs_addr0(i) == calAddrs_addr2(ireuse, jreuse) ) {
                    if ( checkIntervens_addr0_2(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCals_addr0_2(i, ireuse, jreuse) );
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
void pairs_addr1_0() {
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
    break;
}
}
void pairs_addr1_1() {
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
                    if ( calAddrs_addr1(i, j) == calAddrs_addr1(ireuse, jreuse) ) {
                        if ( checkIntervens_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCals_addr1_1(i, j, ireuse, jreuse) );
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
void pairs_addr1_2() {
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
                    if ( calAddrs_addr1(i, j) == calAddrs_addr2(ireuse, jreuse) ) {
                        if ( checkIntervens_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCals_addr1_2(i, j, ireuse, jreuse) );
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
void pairs_addr2_0() {
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
    break;
}
}
void pairs_addr2_1() {
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
                    if ( calAddrs_addr2(i, j) == calAddrs_addr1(ireuse, jreuse) ) {
                        if ( checkIntervens_addr2_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCals_addr2_1(i, j, ireuse, jreuse) );
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
void pairs_addr2_2() {
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
                    if ( calAddrs_addr2(i, j) == calAddrs_addr2(ireuse, jreuse) ) {
                        if ( checkIntervens_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCals_addr2_2(i, j, ireuse, jreuse) );
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
    pairp_addr0_0();
    pairq_addr0_0();
    pairq_addr0_1();
    pairq_addr0_2();
    pairq_addr1_0();
    pairq_addr1_1();
    pairq_addr1_2();
    pairq_addr2_0();
    pairq_addr2_1();
    pairq_addr2_2();
    pairr_addr0_0();
    pairs_addr0_0();
    pairs_addr0_1();
    pairs_addr0_2();
    pairs_addr1_0();
    pairs_addr1_1();
    pairs_addr1_2();
    pairs_addr2_0();
    pairs_addr2_1();
    pairs_addr2_2();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
bicg_cpu */ 
