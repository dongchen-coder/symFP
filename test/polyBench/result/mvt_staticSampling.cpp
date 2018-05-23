
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
------array access x1.addr i
------array access a.addr ((i * 1024) + j)
------array access y1.addr j
------array access x1.addr i
--i
--Loop Bound: (0, 1024)
----j
----Loop Bound: (0, 1024)
------array access x2.addr i
------array access a.addr ((j * 1024) + i)
------array access y2.addr j
------array access x2.addr i

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
int calAddrx1_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddra_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddry1_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrx1_addr1( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrx2_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddra_addr1( int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 32;
    return result;
}
int calAddry2_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrx2_addr1( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
void ref_a_addr0() {
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
                    if ( calAddra_addr0( i, j) == calAddra_addr0(i_Start, j_Start)) {
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
void ref_a_addr1() {
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
                    if ( calAddra_addr1( i, j) == calAddra_addr1(i_Start, j_Start)) {
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
void ref_x1_addr0() {
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
                    if ( calAddrx1_addr0( i, j) == calAddrx1_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx1_addr1( i, j) == calAddrx1_addr0(i_Start, j_Start)) {
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
void ref_x1_addr1() {
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
                    if ( calAddrx1_addr0( i, j) == calAddrx1_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx1_addr1( i, j) == calAddrx1_addr1(i_Start, j_Start)) {
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
void ref_x2_addr0() {
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
                    if ( calAddrx2_addr0( i, j) == calAddrx2_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx2_addr1( i, j) == calAddrx2_addr0(i_Start, j_Start)) {
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
void ref_x2_addr1() {
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
                    if ( calAddrx2_addr0( i, j) == calAddrx2_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx2_addr1( i, j) == calAddrx2_addr1(i_Start, j_Start)) {
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
void ref_y1_addr0() {
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
                    if ( calAddry1_addr0( i, j) == calAddry1_addr0(i_Start, j_Start)) {
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
void ref_y2_addr0() {
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
                    if ( calAddry2_addr0( i, j) == calAddry2_addr0(i_Start, j_Start)) {
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
int main() {
    ref_a_addr0();
    ref_a_addr1();
    ref_x1_addr0();
    ref_x1_addr1();
    ref_x2_addr0();
    ref_x2_addr1();
    ref_y1_addr0();
    ref_y2_addr0();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
runMvt */ 
