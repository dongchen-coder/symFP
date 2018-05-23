
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
----array access s.addr i
--i
--Loop Bound: (0, 1024)
----array access q.addr i
----j
----Loop Bound: (0, 1024)
------array access s.addr j
------array access r.addr i
------array access A.addr ((i * 1024) + j)
------array access s.addr j
------array access q.addr i
------array access A.addr ((i * 1024) + j)
------array access p.addr j
------array access q.addr i

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
int calAddrs_addr0( int i) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrq_addr0( int i) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrs_addr1( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrr_addr0( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrA_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddrs_addr2( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrq_addr1( int i, int j) {
    int result = (i) * 8 / 32;
    return result;
}
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 32;
    return result;
}
int calAddrp_addr0( int i, int j) {
    int result = (j) * 8 / 32;
    return result;
}
int calAddrq_addr2( int i, int j) {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
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
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( i, j) == calAddrA_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( i, j) == calAddrA_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
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
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
            }
        }
EndSample:
        s++;
        }
}
void ref_p_addr0() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( i, j) == calAddrp_addr0(i_Start, j_Start)) {
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
void ref_q_addr0() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr0(i_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            cntStart = true;
            int jLB2 = 0;
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr0(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr0(i_Start)) {
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
void ref_q_addr1() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr1(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr1(i_Start, j_Start)) {
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
void ref_q_addr2() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrq_addr0( i) == calAddrq_addr2(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr1( i, j) == calAddrq_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr2( i, j) == calAddrq_addr2(i_Start, j_Start)) {
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
void ref_r_addr0() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrr_addr0( i, j) == calAddrr_addr0(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
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
void ref_s_addr0() {
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
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            if (cntStart == true) {
                cnt++;
                if ( calAddrs_addr0( i) == calAddrs_addr0(i_Start)) {
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
void ref_s_addr1() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr1( i, j) == calAddrs_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr2( i, j) == calAddrs_addr1(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
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
void ref_s_addr2() {
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
        int iLB1 = i_Start;
        for ( int i = iLB1; i < 1024; i++) {
            if (cntStart == true) cnt++;
            int jLB2 = 0;
            if ( i == i_Start ) {
                jLB2 = j_Start;
            }
            for ( int j = jLB2; j < 1024; j++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr1( i, j) == calAddrs_addr2(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrs_addr2( i, j) == calAddrs_addr2(i_Start, j_Start)) {
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_p_addr0();
    ref_q_addr0();
    ref_q_addr1();
    ref_q_addr2();
    ref_r_addr0();
    ref_s_addr0();
    ref_s_addr1();
    ref_s_addr2();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
bicg_cpu */ 
