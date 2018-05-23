
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
------array access A.addr ((i * 1024) + j)
------array access u1.addr i
------array access v1.addr j
------array access u2.addr i
------array access v2.addr j
------array access A.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access x.addr i
------array access A.addr ((j * 1024) + i)
------array access y.addr j
------array access x.addr i
--i
--Loop Bound: (0, 1024)
----array access x.addr i
----array access z.addr i
----array access x.addr i
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access w.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access w.addr i

Finish analysis loops */ 
 // Start to generating Static Sampling Code (reference based)
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
int calAddrA_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddru1_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrv1_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddru2_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrv2_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddrx_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrA_addr2( int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 32;
    return result;
}
int calAddry_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrx_addr1( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrx_addr2( int i) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrz_addr0( int i) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrx_addr3( int i) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrw_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrA_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddrx_addr4( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrw_addr1( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
void ref_A_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
        }
EndSample:
        s++;
        }
}
void ref_A_addr1() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
        }
EndSample:
        s++;
        }
}
void ref_A_addr2() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( i, j) == calAddrA_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_A_addr3() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( i, j) == calAddrA_addr3(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_u1_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru1_addr0( i, j) == calAddru1_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_u2_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru2_addr0( i, j) == calAddru2_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_v1_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv1_addr0( i, j) == calAddrv1_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_v2_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv2_addr0( i, j) == calAddrv2_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_w_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( i, j) == calAddrw_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( i, j) == calAddrw_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
        }
EndSample:
        s++;
        }
}
void ref_w_addr1() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( i, j) == calAddrw_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( i, j) == calAddrw_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
        }
EndSample:
        s++;
        }
}
void ref_x_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( i, j) == calAddrx_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
        }
EndSample:
        s++;
        }
}
void ref_x_addr1() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( i, j) == calAddrx_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
        }
EndSample:
        s++;
        }
}
void ref_x_addr2() {
    set<string> record;
    for ( int s = 0; s < 20;) {
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr2(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr2(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
        }
EndSample:
        s++;
        }
}
void ref_x_addr3() {
    set<string> record;
    for ( int s = 0; s < 20;) {
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( i) == calAddrx_addr3(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( i) == calAddrx_addr3(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
        }
EndSample:
        s++;
        }
}
void ref_x_addr4() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB5 = i_Start;
        for ( int i = iLB5; i < 1024; i++) {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( i, j) == calAddrx_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_y_addr0() {
    set<string> record;
    for ( int s = 0; s < 400;) {
        int i_Start = rand() % (1024 - 0) + 0;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            j_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB2 = i_Start;
        for ( int i = iLB2; i < 1024; i++) {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( i, j) == calAddry_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_z_addr0() {
    set<string> record;
    for ( int s = 0; s < 20;) {
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        while ( record.find(idx_string) != record.end() ) {
            i_Start = rand() % (1024 - 0) + 0;
            idx_string = std::to_string(i_Start) + "_" ;
        }
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;
        int iLB4 = i_Start;
        for ( int i = iLB4; i < 1024; i++) {
            if (cntStart == true) cnt++;
            if (cntStart == true) {
                cnt++;
                if ( calAddrz_addr0( i) == calAddrz_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) cnt++;
        }
EndSample:
        s++;
        }
}
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr3();
    ref_u1_addr0();
    ref_u2_addr0();
    ref_v1_addr0();
    ref_v2_addr0();
    ref_w_addr0();
    ref_w_addr1();
    ref_x_addr0();
    ref_x_addr1();
    ref_x_addr2();
    ref_x_addr3();
    ref_x_addr4();
    ref_y_addr0();
    ref_z_addr0();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
gemver */ 
