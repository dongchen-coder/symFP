
 /* Start to analysis array index
Array index info
a.addr ((1026 * i) + j)
a.addr (((1026 * i) + j) + 1)
a.addr (((1026 * i) + j) - 1)
a.addr ((1026 * (i - 1)) + j)
a.addr ((1026 * (i + 1)) + j)
b.addr ((1026 * i) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %a
double* %b

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (1, 1025)
----j
----Loop Bound: (1, 1025)
------array access   %5 = load double, double* %arrayidx, align 8
------array access   %9 = load double, double* %arrayidx8, align 8
------array access   %13 = load double, double* %arrayidx13, align 8
------array access   %17 = load double, double* %arrayidx19, align 8
------array access   %21 = load double, double* %arrayidx25, align 8
------array access   store double %add26, double* %arrayidx30, align 8

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
int calAddra_addr0( int i, int j) {
    int result = ((1026 * i) + j);
    return result;
}
int calAddra_addr1( int i, int j) {
    int result = (((1026 * i) + j) + 1);
    return result;
}
int calAddra_addr2( int i, int j) {
    int result = (((1026 * i) + j) - 1);
    return result;
}
int calAddra_addr3( int i, int j) {
    int result = ((1026 * (i - 1)) + j);
    return result;
}
int calAddra_addr4( int i, int j) {
    int result = ((1026 * (i + 1)) + j);
    return result;
}
int calAddrb_addr0( int i, int j) {
    int result = ((1026 * i) + j);
    return result;
}
int rtCala_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 0;
}
int rtCala_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 0;
}
int rtCala_addr0_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 0;
}
int rtCala_addr0_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 0;
}
int rtCala_addr0_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 0;
}
int rtCala_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 1;
}
int rtCala_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 1;
}
int rtCala_addr1_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 1;
}
int rtCala_addr1_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 1;
}
int rtCala_addr1_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 1;
}
int rtCala_addr2_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 2;
}
int rtCala_addr2_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 2;
}
int rtCala_addr2_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 2;
}
int rtCala_addr2_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 2;
}
int rtCala_addr2_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 2;
}
int rtCala_addr3_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 3;
}
int rtCala_addr3_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 3;
}
int rtCala_addr3_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 3;
}
int rtCala_addr3_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 3;
}
int rtCala_addr3_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 3;
}
int rtCala_addr4_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 0 - 4;
}
int rtCala_addr4_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 1 - 4;
}
int rtCala_addr4_2(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 2 - 4;
}
int rtCala_addr4_3(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 3 - 4;
}
int rtCala_addr4_4(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 4 - 4;
}
int rtCalb_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 6144 + (jreuse - j) * 6 + 5 - 5;
}
bool checkIntervena_addr0_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr0_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr0_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr0_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr0_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr0(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr1(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr2_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr2_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr2_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr2_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr2_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr2(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr2(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr3_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr3_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr3_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr3_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr3_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr3(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr3(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr4_0(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr4(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr4_1(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr4(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr4_2(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr4(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr4_3(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse - 1 ;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr4(iInterven, jInterven) == calAddra_addr4(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr4_4(int i, int j, int ireuse, int jreuse) { 
    for(int iInterven = i; iInterven <= ireuse; iInterven++) {
        int jInterven;
        int jIntervenUB;
        if (iInterven == i) {
            jInterven = j + 1 ;
        } else {
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr1(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr2(iInterven, jInterven) == calAddra_addr4(i, j)) {
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
            jInterven = 1;
        }
        if (iInterven == ireuse) {
            jIntervenUB = jreuse;
        } else {
            jIntervenUB = 1025- 1;
        }
        for(; jInterven <= jIntervenUB; jInterven++) {
            if( calAddra_addr3(iInterven, jInterven) == calAddra_addr4(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenb_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
void paira_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr0(i, j) == calAddra_addr0(ireuse, jreuse) ) {
                        if ( checkIntervena_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr0_0(i, j, ireuse, jreuse) );
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
void paira_addr0_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr0(i, j) == calAddra_addr1(ireuse, jreuse) ) {
                        if ( checkIntervena_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr0_1(i, j, ireuse, jreuse) );
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
void paira_addr0_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr0(i, j) == calAddra_addr2(ireuse, jreuse) ) {
                        if ( checkIntervena_addr0_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr0_2(i, j, ireuse, jreuse) );
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
void paira_addr0_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr0(i, j) == calAddra_addr3(ireuse, jreuse) ) {
                        if ( checkIntervena_addr0_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr0_3(i, j, ireuse, jreuse) );
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
void paira_addr0_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr0(i, j) == calAddra_addr4(ireuse, jreuse) ) {
                        if ( checkIntervena_addr0_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr0_4(i, j, ireuse, jreuse) );
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
void paira_addr1_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr1(i, j) == calAddra_addr0(ireuse, jreuse) ) {
                        if ( checkIntervena_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr1_0(i, j, ireuse, jreuse) );
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
void paira_addr1_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr1(i, j) == calAddra_addr1(ireuse, jreuse) ) {
                        if ( checkIntervena_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr1_1(i, j, ireuse, jreuse) );
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
void paira_addr1_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr1(i, j) == calAddra_addr2(ireuse, jreuse) ) {
                        if ( checkIntervena_addr1_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr1_2(i, j, ireuse, jreuse) );
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
void paira_addr1_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr1(i, j) == calAddra_addr3(ireuse, jreuse) ) {
                        if ( checkIntervena_addr1_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr1_3(i, j, ireuse, jreuse) );
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
void paira_addr1_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr1(i, j) == calAddra_addr4(ireuse, jreuse) ) {
                        if ( checkIntervena_addr1_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr1_4(i, j, ireuse, jreuse) );
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
void paira_addr2_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr2(i, j) == calAddra_addr0(ireuse, jreuse) ) {
                        if ( checkIntervena_addr2_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr2_0(i, j, ireuse, jreuse) );
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
void paira_addr2_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr2(i, j) == calAddra_addr1(ireuse, jreuse) ) {
                        if ( checkIntervena_addr2_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr2_1(i, j, ireuse, jreuse) );
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
void paira_addr2_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr2(i, j) == calAddra_addr2(ireuse, jreuse) ) {
                        if ( checkIntervena_addr2_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr2_2(i, j, ireuse, jreuse) );
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
void paira_addr2_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr2(i, j) == calAddra_addr3(ireuse, jreuse) ) {
                        if ( checkIntervena_addr2_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr2_3(i, j, ireuse, jreuse) );
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
void paira_addr2_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr2(i, j) == calAddra_addr4(ireuse, jreuse) ) {
                        if ( checkIntervena_addr2_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr2_4(i, j, ireuse, jreuse) );
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
void paira_addr3_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr3(i, j) == calAddra_addr0(ireuse, jreuse) ) {
                        if ( checkIntervena_addr3_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr3_0(i, j, ireuse, jreuse) );
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
void paira_addr3_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr3(i, j) == calAddra_addr1(ireuse, jreuse) ) {
                        if ( checkIntervena_addr3_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr3_1(i, j, ireuse, jreuse) );
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
void paira_addr3_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr3(i, j) == calAddra_addr2(ireuse, jreuse) ) {
                        if ( checkIntervena_addr3_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr3_2(i, j, ireuse, jreuse) );
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
void paira_addr3_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr3(i, j) == calAddra_addr3(ireuse, jreuse) ) {
                        if ( checkIntervena_addr3_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr3_3(i, j, ireuse, jreuse) );
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
void paira_addr3_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr3(i, j) == calAddra_addr4(ireuse, jreuse) ) {
                        if ( checkIntervena_addr3_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr3_4(i, j, ireuse, jreuse) );
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
void paira_addr4_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr4(i, j) == calAddra_addr0(ireuse, jreuse) ) {
                        if ( checkIntervena_addr4_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr4_0(i, j, ireuse, jreuse) );
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
void paira_addr4_1() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr4(i, j) == calAddra_addr1(ireuse, jreuse) ) {
                        if ( checkIntervena_addr4_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr4_1(i, j, ireuse, jreuse) );
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
void paira_addr4_2() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr4(i, j) == calAddra_addr2(ireuse, jreuse) ) {
                        if ( checkIntervena_addr4_2(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr4_2(i, j, ireuse, jreuse) );
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
void paira_addr4_3() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr4(i, j) == calAddra_addr3(ireuse, jreuse) ) {
                        if ( checkIntervena_addr4_3(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr4_3(i, j, ireuse, jreuse) );
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
void paira_addr4_4() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddra_addr4(i, j) == calAddra_addr4(ireuse, jreuse) ) {
                        if ( checkIntervena_addr4_4(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCala_addr4_4(i, j, ireuse, jreuse) );
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
void pairb_addr0_0() {
    set<string> record;
    for ( int s = 0; s < 1000; s++) {
        int i = rand() % (1025 - 1) + 1;
        int j = rand() % (1025 - 1) + 1;
        string idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i = rand() % (1025 - 1) + 1;
            j = rand() % (1025 - 1) + 1;
            idx_string = std::to_string(i) + "_" + std::to_string(j) + "_" ;
        }
        record.insert( idx_string );
        bool findReuseFlag = false;
        for ( int ireuse = i; ireuse < 1025; ireuse++) {
            int jreuse;
            if (ireuse == i) {
                jreuse = j + 1;
            } else {
                jreuse = 1;
            }
            for ( ; jreuse < 1025; jreuse++) {
                    if ( calAddrb_addr0(i, j) == calAddrb_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenb_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalb_addr0_0(i, j, ireuse, jreuse) );
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
    paira_addr0_0();
    paira_addr0_1();
    paira_addr0_2();
    paira_addr0_3();
    paira_addr0_4();
    paira_addr1_0();
    paira_addr1_1();
    paira_addr1_2();
    paira_addr1_3();
    paira_addr1_4();
    paira_addr2_0();
    paira_addr2_1();
    paira_addr2_2();
    paira_addr2_3();
    paira_addr2_4();
    paira_addr3_0();
    paira_addr3_1();
    paira_addr3_2();
    paira_addr3_3();
    paira_addr3_4();
    paira_addr4_0();
    paira_addr4_1();
    paira_addr4_2();
    paira_addr4_3();
    paira_addr4_4();
    pairb_addr0_0();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
stencil */ 
