
 /* Start to analysis array index
Array index info
A.addr ((i * 1024) + j)
u1.addr i
v1.addr j
u2.addr i
v2.addr j
A.addr ((i * 1024) + j)
x.addr i
A.addr ((j * 1024) + i)
y.addr j
x.addr i
x.addr i
z.addr i
x.addr i
w.addr i
A.addr ((i * 1024) + j)
x.addr j
w.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double %alpha
double %beta
double* %A
double* %u1
double* %v1
double* %u2
double* %v2
double* %w
double* %x
double* %y
double* %z

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %5 = load double, double* %arrayidx, align 8
------array access   %8 = load double, double* %arrayidx5, align 8
------array access   %11 = load double, double* %arrayidx7, align 8
------array access   %14 = load double, double* %arrayidx11, align 8
------array access   %17 = load double, double* %arrayidx13, align 8
------array access   store double %add15, double* %arrayidx19, align 8
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %27 = load double, double* %arrayidx30, align 8
------array access   %32 = load double, double* %arrayidx34, align 8
------array access   %35 = load double, double* %arrayidx37, align 8
------array access   store double %add39, double* %arrayidx41, align 8
--i
--Loop Bound: (0, 1024)
----array access   %43 = load double, double* %arrayidx52, align 8
----array access   %46 = load double, double* %arrayidx54, align 8
----array access   store double %add55, double* %arrayidx57, align 8
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %54 = load double, double* %arrayidx68, align 8
------array access   %59 = load double, double* %arrayidx72, align 8
------array access   %62 = load double, double* %arrayidx75, align 8
------array access   store double %add77, double* %arrayidx79, align 8

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
int calAddrA_addr0( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddru1_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddrv1_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddru2_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddrv2_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrx_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddrA_addr2( int i, int j) {
    int result = ((j * 1024) + i);
    return result;
}
int calAddry_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrx_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddrx_addr2( int i) {
    int result = i;
    return result;
}
int calAddrz_addr0( int i) {
    int result = i;
    return result;
}
int calAddrx_addr3( int i) {
    int result = i;
    return result;
}
int calAddrw_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddrA_addr3( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddrx_addr4( int i, int j) {
    int result = j;
    return result;
}
int calAddrw_addr1( int i, int j) {
    int result = i;
    return result;
}
int rtCalA_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 0;
}
int rtCalA_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 5 - 0;
}
int rtCalA_addr0_2(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 6144 + (ireuse - 0) * 4096 + (1024 - j) * 6 + (jreuse - 0) * 4 + 1 - 0;
}
int rtCalA_addr0_3(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 6144 + (ireuse - 0) * 4096 + (1024 - j) * 6 + (jreuse - 0) * 4 + 4096 + 3 + 1 - 0;
}
int rtCalA_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 5;
}
int rtCalA_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 5 - 5;
}
int rtCalA_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 6144 + (ireuse - 0) * 4096 + (1024 - j) * 6 + (jreuse - 0) * 4 + 1 - 5;
}
int rtCalA_addr1_3(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 6144 + (ireuse - 0) * 4096 + (1024 - j) * 6 + (jreuse - 0) * 4 + 4096 + 3 + 1 - 5;
}
int rtCalA_addr2_0(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 6144 + (1024 - j) * 4 + (jreuse - 0) * 6 + 0 - 1;
}
int rtCalA_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 6144 + (1024 - j) * 4 + (jreuse - 0) * 6 + 5 - 1;
}
int rtCalA_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 1 - 1;
}
int rtCalA_addr2_3(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 3 + 1 - 1;
}
int rtCalA_addr3_0(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 6144 + (1024 - j) * 4 + (jreuse - 0) * 6 + 0 - 1;
}
int rtCalA_addr3_1(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 6144 + (1024 - j) * 4 + (jreuse - 0) * 6 + 5 - 1;
}
int rtCalA_addr3_2(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 1 - 1;
}
int rtCalA_addr3_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 1 - 1;
}
int rtCalu1_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 1;
}
int rtCalu2_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 3;
}
int rtCalv1_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 2;
}
int rtCalv2_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 4;
}
int rtCalw_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 0;
}
int rtCalw_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 0;
}
int rtCalw_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 3;
}
int rtCalw_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 3;
}
int rtCalx_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 0;
}
int rtCalx_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 0;
}
int rtCalx_addr0_2(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 0 - 0;
}
int rtCalx_addr0_3(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 2 - 0;
}
int rtCalx_addr0_4(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 3 + 2 - 0;
}
int rtCalx_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 3;
}
int rtCalx_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 3;
}
int rtCalx_addr1_2(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 0 - 3;
}
int rtCalx_addr1_3(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 2 - 3;
}
int rtCalx_addr1_4(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 3 + 2 - 3;
}
int rtCalx_addr2_0(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 0 - 0;
}
int rtCalx_addr2_1(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 3 - 0;
}
int rtCalx_addr2_2(int i, int ireuse) {
    return (ireuse - i) * 3 + 0 - 0;
}
int rtCalx_addr2_3(int i, int ireuse) {
    return (ireuse - i) * 3 + 2 - 0;
}
int rtCalx_addr2_4(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 2 - 0;
}
int rtCalx_addr3_0(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 0 - 2;
}
int rtCalx_addr3_1(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 3 - 2;
}
int rtCalx_addr3_2(int i, int ireuse) {
    return (ireuse - i) * 3 + 0 - 2;
}
int rtCalx_addr3_3(int i, int ireuse) {
    return (ireuse - i) * 3 + 2 - 2;
}
int rtCalx_addr3_4(int i, int ireuse, int jreuse) {
    return (1024 - i) * 3 + (ireuse - 0) * 4096 + (jreuse - 0) * 4 + 2 - 2;
}
int rtCalx_addr4_0(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 0 - 2;
}
int rtCalx_addr4_1(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 3 - 2;
}
int rtCalx_addr4_2(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 0 - 2;
}
int rtCalx_addr4_3(int i, int j, int ireuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 3 + (1024 - j) * 4 + 2 - 2;
}
int rtCalx_addr4_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 2 - 2;
}
int rtCaly_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 2 - 2;
}
int rtCalz_addr0_0(int i, int ireuse) {
    return (ireuse - i) * 3 + 1 - 1;
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
bool checkIntervenA_addr0_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
bool checkIntervenA_addr0_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr2(iInterven, jInterven) == calAddrA_addr0(i, j)) {
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
bool checkIntervenA_addr1_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
bool checkIntervenA_addr1_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
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
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == 0) {
            jInterven = 0;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr2(iInterven, jInterven) == calAddrA_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr2_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr2_1(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr2_2(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr2_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrA_addr2(iInterven, jInterven) == calAddrA_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenA_addr3_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr3_1(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr3_2(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenA_addr3_3(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenu1_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenu2_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenv1_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenv2_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenw_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrw_addr1(iInterven, jInterven) == calAddrw_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenw_addr0_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrw_addr0(iInterven, jInterven) == calAddrw_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenw_addr1_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrw_addr1(iInterven, jInterven) == calAddrw_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenw_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrw_addr0(iInterven, jInterven) == calAddrw_addr1(i, j)) {
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
bool checkIntervenx_addr0_2(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr0(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr0_3(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr0(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr0_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr0(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr0(i, j)) {
            return true;
        }
    }
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr0(i, j)) {
            return true;
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
bool checkIntervenx_addr1_2(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr1_3(int i, int j, int ireuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= ireuse; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr1_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr0(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 0;
        }
        if (iInterven == 1024) {
            jIntervenUB = 1024;
        } else {
            jIntervenUB = 1024- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddrx_addr1(iInterven, jInterven) == calAddrx_addr1(i, j)) {
                return true;
            }
        }
    }
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr1(i, j)) {
            return true;
        }
    }
    for(int iInterven = 0; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr1(i, j)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr2_0(int i, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr2_1(int i, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr2_2(int i, int ireuse) { 
    for(int iInterven = i; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr2(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr2_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr2(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr2_4(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr2(i)) {
            return true;
        }
    }
    for(int iInterven = i; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr2(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr3_0(int i, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr3_1(int i, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr3_2(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse - 1 ; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr3_3(int i, int ireuse) { 
    for(int iInterven = i + 1 ; iInterven <= ireuse; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr3_4(int i, int ireuse, int jreuse) { 
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr2(iInterven) == calAddrx_addr3(i)) {
            return true;
        }
    }
    for(int iInterven = i + 1 ; iInterven <= 1024; iInterven++) {
        if( calAddrx_addr3(iInterven) == calAddrx_addr3(i)) {
            return true;
        }
    }
    return false;
}
bool checkIntervenx_addr4_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr4_1(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx_addr4_2(int i, int j, int ireuse) { 
    return false;
}
bool checkIntervenx_addr4_3(int i, int j, int ireuse) { 
    return false;
}
bool checkIntervenx_addr4_4(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkInterveny_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenz_addr0_0(int i, int ireuse) { 
    return false;
}
void pairA_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr0_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr0(i, j) == calAddrA_addr2(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr0_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr0_2(i, j, ireuse, jreuse) );
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
void pairA_addr0_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr0(i, j) == calAddrA_addr3(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr0_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr0_3(i, j, ireuse, jreuse) );
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
    for ( int s = 0; s < 1000; s++) {
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
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr1(i, j) == calAddrA_addr2(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr1_2(i, j, ireuse, jreuse) );
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
void pairA_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr1(i, j) == calAddrA_addr3(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr1_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr1_3(i, j, ireuse, jreuse) );
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
void pairA_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrA_addr2(i, j) == calAddrA_addr2(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr2_2(i, j, ireuse, jreuse) );
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
void pairA_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrA_addr2(i, j) == calAddrA_addr3(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr2_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr2_3(i, j, ireuse, jreuse) );
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
void pairA_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairA_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrA_addr3(i, j) == calAddrA_addr3(ireuse, jreuse) ) {
                        if ( checkIntervenA_addr3_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalA_addr3_3(i, j, ireuse, jreuse) );
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
void pairu1_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddru1_addr0(i, j) == calAddru1_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenu1_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalu1_addr0_0(i, j, ireuse, jreuse) );
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
void pairu2_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddru2_addr0(i, j) == calAddru2_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenu2_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalu2_addr0_0(i, j, ireuse, jreuse) );
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
void pairv1_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrv1_addr0(i, j) == calAddrv1_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenv1_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalv1_addr0_0(i, j, ireuse, jreuse) );
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
void pairv2_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrv2_addr0(i, j) == calAddrv2_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenv2_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalv2_addr0_0(i, j, ireuse, jreuse) );
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
void pairw_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrw_addr0(i, j) == calAddrw_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenw_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalw_addr0_0(i, j, ireuse, jreuse) );
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
void pairw_addr0_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrw_addr0(i, j) == calAddrw_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenw_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalw_addr0_1(i, j, ireuse, jreuse) );
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
void pairw_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrw_addr1(i, j) == calAddrw_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenw_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalw_addr1_0(i, j, ireuse, jreuse) );
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
void pairw_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrw_addr1(i, j) == calAddrw_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenw_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalw_addr1_1(i, j, ireuse, jreuse) );
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
void pairx_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr0_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
                if ( calAddrx_addr0(i, j) == calAddrx_addr2(ireuse) ) {
                    if ( checkIntervenx_addr0_2(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalx_addr0_2(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairx_addr0_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
                if ( calAddrx_addr0(i, j) == calAddrx_addr3(ireuse) ) {
                    if ( checkIntervenx_addr0_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalx_addr0_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairx_addr0_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr0(i, j) == calAddrx_addr4(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr0_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr0_4(i, j, ireuse, jreuse) );
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
    for ( int s = 0; s < 1000; s++) {
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
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
                if ( calAddrx_addr1(i, j) == calAddrx_addr2(ireuse) ) {
                    if ( checkIntervenx_addr1_2(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalx_addr1_2(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairx_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
                if ( calAddrx_addr1(i, j) == calAddrx_addr3(ireuse) ) {
                    if ( checkIntervenx_addr1_3(i, j, ireuse) == false) {
                      rtHistoCal(  rtCalx_addr1_3(i, j, ireuse) );
                    }
                    findReuseFlag = true;
                }
                if (findReuseFlag == true) {
                    break;
                }
            }
        }
}
void pairx_addr1_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
        for ( int ireuse = 0; ireuse < 1024; ireuse++) {
            for ( int jreuse = 0; jreuse < 1024; jreuse++) {
                    if ( calAddrx_addr1(i, j) == calAddrx_addr4(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr1_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr1_4(i, j, ireuse, jreuse) );
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
void pairx_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairx_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairx_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrx_addr2(i) == calAddrx_addr2(ireuse) ) {
                if ( checkIntervenx_addr2_2(i, ireuse) == false) {
                  rtHistoCal(  rtCalx_addr2_2(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairx_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1024; ireuse++) {
            if ( calAddrx_addr2(i) == calAddrx_addr3(ireuse) ) {
                if ( checkIntervenx_addr2_3(i, ireuse) == false) {
                  rtHistoCal(  rtCalx_addr2_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairx_addr2_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                if ( calAddrx_addr2(i) == calAddrx_addr4(ireuse, jreuse) ) {
                    if ( checkIntervenx_addr2_4(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCalx_addr2_4(i, ireuse, jreuse) );
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
void pairx_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairx_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
    break;
}
}
void pairx_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrx_addr3(i) == calAddrx_addr2(ireuse) ) {
                if ( checkIntervenx_addr3_2(i, ireuse) == false) {
                  rtHistoCal(  rtCalx_addr3_2(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairx_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrx_addr3(i) == calAddrx_addr3(ireuse) ) {
                if ( checkIntervenx_addr3_3(i, ireuse) == false) {
                  rtHistoCal(  rtCalx_addr3_3(i, ireuse) );
                }
                findReuseFlag = true;
            }
            if (findReuseFlag == true) {
                break;
            }
        }
    }
}
void pairx_addr3_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                if ( calAddrx_addr3(i) == calAddrx_addr4(ireuse, jreuse) ) {
                    if ( checkIntervenx_addr3_4(i, ireuse, jreuse) == false) {
                      rtHistoCal(  rtCalx_addr3_4(i, ireuse, jreuse) );
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
void pairx_addr4_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr4_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr4_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr4_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
void pairx_addr4_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddrx_addr4(i, j) == calAddrx_addr4(ireuse, jreuse) ) {
                        if ( checkIntervenx_addr4_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx_addr4_4(i, j, ireuse, jreuse) );
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
    for ( int s = 0; s < 1000; s++) {
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
                    if ( calAddry_addr0(i, j) == calAddry_addr0(ireuse, jreuse) ) {
                        if ( checkInterveny_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly_addr0_0(i, j, ireuse, jreuse) );
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
void pairz_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i + 1; ireuse < 1024; ireuse++) {
            if ( calAddrz_addr0(i) == calAddrz_addr0(ireuse) ) {
                if ( checkIntervenz_addr0_0(i, ireuse) == false) {
                  rtHistoCal(  rtCalz_addr0_0(i, ireuse) );
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
    pairA_addr0_1();
    pairA_addr0_2();
    pairA_addr0_3();
    pairA_addr1_0();
    pairA_addr1_1();
    pairA_addr1_2();
    pairA_addr1_3();
    pairA_addr2_0();
    pairA_addr2_1();
    pairA_addr2_2();
    pairA_addr2_3();
    pairA_addr3_0();
    pairA_addr3_1();
    pairA_addr3_2();
    pairA_addr3_3();
    pairu1_addr0_0();
    pairu2_addr0_0();
    pairv1_addr0_0();
    pairv2_addr0_0();
    pairw_addr0_0();
    pairw_addr0_1();
    pairw_addr1_0();
    pairw_addr1_1();
    pairx_addr0_0();
    pairx_addr0_1();
    pairx_addr0_2();
    pairx_addr0_3();
    pairx_addr0_4();
    pairx_addr1_0();
    pairx_addr1_1();
    pairx_addr1_2();
    pairx_addr1_3();
    pairx_addr1_4();
    pairx_addr2_0();
    pairx_addr2_1();
    pairx_addr2_2();
    pairx_addr2_3();
    pairx_addr2_4();
    pairx_addr3_0();
    pairx_addr3_1();
    pairx_addr3_2();
    pairx_addr3_3();
    pairx_addr3_4();
    pairx_addr4_0();
    pairx_addr4_1();
    pairx_addr4_2();
    pairx_addr4_3();
    pairx_addr4_4();
    pairy_addr0_0();
    pairz_addr0_0();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
gemver */ 
