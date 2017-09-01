
 /* Start to analysis array index
Array index info
a.addr ((i * 1024) + j)
x1.addr i
y1.addr j
x1.addr i
x2.addr i
a.addr ((j * 1024) + i)
y2.addr j
x2.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %n
double* %a
double* %x1
double* %x2
double* %y1
double* %y2

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %4 = load double, double* %arrayidx, align 8
------array access   %8 = load double, double* %arrayidx5, align 8
------array access   %11 = load double, double* %arrayidx7, align 8
------array access   store double %add9, double* %arrayidx11, align 8
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access   %20 = load double, double* %arrayidx22, align 8
------array access   %24 = load double, double* %arrayidx26, align 8
------array access   %27 = load double, double* %arrayidx28, align 8
------array access   store double %add30, double* %arrayidx32, align 8

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
int calAddrx1_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddra_addr0( int i, int j) {
    int result = ((i * 1024) + j);
    return result;
}
int calAddry1_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrx1_addr1( int i, int j) {
    int result = i;
    return result;
}
int calAddrx2_addr0( int i, int j) {
    int result = i;
    return result;
}
int calAddra_addr1( int i, int j) {
    int result = ((j * 1024) + i);
    return result;
}
int calAddry2_addr0( int i, int j) {
    int result = j;
    return result;
}
int calAddrx2_addr1( int i, int j) {
    int result = i;
    return result;
}
int rtCala_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 1 - 1;
}
int rtCala_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 1 - 1;
}
int rtCala_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (1024 - i) * 4096 + (ireuse - 0) * 4096 + (1024 - j) * 4 + (jreuse - 0) * 4 + 1 - 1;
}
int rtCala_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 1 - 1;
}
int rtCalx1_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 0;
}
int rtCalx1_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 0;
}
int rtCalx1_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 3;
}
int rtCalx1_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 3;
}
int rtCalx2_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 0;
}
int rtCalx2_addr0_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 0;
}
int rtCalx2_addr1_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 0 - 3;
}
int rtCalx2_addr1_1(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 3 - 3;
}
int rtCaly1_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 2 - 2;
}
int rtCaly2_addr0_0(int i, int j, int ireuse, int jreuse) {
    return (ireuse - i) * 4096 + (jreuse - j) * 4 + 2 - 2;
}
bool checkIntervena_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervena_addr0_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddra_addr0(iInterven, jInterven) == calAddra_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervena_addr1_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervena_addr1_1(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkIntervenx1_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx1_addr1(iInterven, jInterven) == calAddrx1_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx1_addr0_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx1_addr0(iInterven, jInterven) == calAddrx1_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx1_addr1_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx1_addr1(iInterven, jInterven) == calAddrx1_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx1_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx1_addr0(iInterven, jInterven) == calAddrx1_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx2_addr0_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx2_addr1(iInterven, jInterven) == calAddrx2_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx2_addr0_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx2_addr0(iInterven, jInterven) == calAddrx2_addr0(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx2_addr1_0(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx2_addr1(iInterven, jInterven) == calAddrx2_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkIntervenx2_addr1_1(int i, int j, int ireuse, int jreuse) { 
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
            if( calAddrx2_addr0(iInterven, jInterven) == calAddrx2_addr1(i, j)) {
                return true;
            }
        }
    }
    return false;
}
bool checkInterveny1_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
bool checkInterveny2_addr0_0(int i, int j, int ireuse, int jreuse) { 
    return false;
}
void paira_addr0_0() {
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
void paira_addr1_0() {
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
void paira_addr1_1() {
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
void pairx1_addr0_0() {
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
                    if ( calAddrx1_addr0(i, j) == calAddrx1_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx1_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx1_addr0_0(i, j, ireuse, jreuse) );
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
void pairx1_addr0_1() {
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
                    if ( calAddrx1_addr0(i, j) == calAddrx1_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx1_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx1_addr0_1(i, j, ireuse, jreuse) );
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
void pairx1_addr1_0() {
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
                    if ( calAddrx1_addr1(i, j) == calAddrx1_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx1_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx1_addr1_0(i, j, ireuse, jreuse) );
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
void pairx1_addr1_1() {
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
                    if ( calAddrx1_addr1(i, j) == calAddrx1_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx1_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx1_addr1_1(i, j, ireuse, jreuse) );
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
void pairx2_addr0_0() {
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
                    if ( calAddrx2_addr0(i, j) == calAddrx2_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx2_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx2_addr0_0(i, j, ireuse, jreuse) );
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
void pairx2_addr0_1() {
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
                    if ( calAddrx2_addr0(i, j) == calAddrx2_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx2_addr0_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx2_addr0_1(i, j, ireuse, jreuse) );
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
void pairx2_addr1_0() {
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
                    if ( calAddrx2_addr1(i, j) == calAddrx2_addr0(ireuse, jreuse) ) {
                        if ( checkIntervenx2_addr1_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx2_addr1_0(i, j, ireuse, jreuse) );
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
void pairx2_addr1_1() {
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
                    if ( calAddrx2_addr1(i, j) == calAddrx2_addr1(ireuse, jreuse) ) {
                        if ( checkIntervenx2_addr1_1(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCalx2_addr1_1(i, j, ireuse, jreuse) );
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
void pairy1_addr0_0() {
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
                    if ( calAddry1_addr0(i, j) == calAddry1_addr0(ireuse, jreuse) ) {
                        if ( checkInterveny1_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly1_addr0_0(i, j, ireuse, jreuse) );
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
void pairy2_addr0_0() {
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
                    if ( calAddry2_addr0(i, j) == calAddry2_addr0(ireuse, jreuse) ) {
                        if ( checkInterveny2_addr0_0(i, j, ireuse, jreuse) == false) {
                          rtHistoCal(  rtCaly2_addr0_0(i, j, ireuse, jreuse) );
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
    paira_addr1_0();
    paira_addr1_1();
    pairx1_addr0_0();
    pairx1_addr0_1();
    pairx1_addr1_0();
    pairx1_addr1_1();
    pairx2_addr0_0();
    pairx2_addr0_1();
    pairx2_addr1_0();
    pairx2_addr1_1();
    pairy1_addr0_0();
    pairy2_addr0_0();
    rtDump();
    RTtoMR_AET();    dumpMR();    return 0;
}
 /* Start to analyze function:  
runMvt */ 
