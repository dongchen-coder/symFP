
 /* Start to analysis array index
Array index info
v.addr (0 + i)
p.addr ((i * 1024) + 0)
v.addr (0 + i)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
u.addr (((j * 1024) + i) - 1)
u.addr ((j * 1024) + i)
u.addr (((j * 1024) + i) + 1)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
v.addr (1047552 + i)
p.addr ((i * 1024) + j)
v.addr (((j + 1) * 1024) + i)
q.addr ((i * 1024) + j)
v.addr ((j * 1024) + i)
u.addr ((i * 1024) + 0)
p.addr ((i * 1024) + 0)
u.addr ((i * 1024) + 0)
q.addr ((i * 1024) + 0)
p.addr (((i * 1024) + j) - 1)
p.addr ((i * 1024) + j)
v.addr (((i - 1) * 1024) + j)
v.addr ((i * 1024) + j)
v.addr (((i + 1) * 1024) + j)
q.addr (((i * 1024) + j) - 1)
p.addr (((i * 1024) + j) - 1)
q.addr ((i * 1024) + j)
u.addr (((i * 1024) + 1024) - 1)
p.addr ((i * 1024) + j)
u.addr (((i * 1024) + j) + 1)
q.addr ((i * 1024) + j)
u.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %p
double* %q
double* %v
double* %u

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access v.addr (0 + i)
------array access p.addr ((i * 1024) + 0)
------array access v.addr (0 + i)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((j * 1024) + i) - 1)
--------array access u.addr ((j * 1024) + i)
--------array access u.addr (((j * 1024) + i) + 1)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access v.addr (1047552 + i)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((j + 1) * 1024) + i)
--------array access q.addr ((i * 1024) + j)
--------array access v.addr ((j * 1024) + i)
----i
----Loop Bound: (1, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------array access u.addr ((i * 1024) + 0)
------array access p.addr ((i * 1024) + 0)
------array access u.addr ((i * 1024) + 0)
------array access q.addr ((i * 1024) + 0)
------j
------Loop Bound: (1, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access p.addr (((i * 1024) + j) - 1)
--------array access p.addr ((i * 1024) + j)
--------array access v.addr (((i - 1) * 1024) + j)
--------array access v.addr ((i * 1024) + j)
--------array access v.addr (((i + 1) * 1024) + j)
--------array access q.addr (((i * 1024) + j) - 1)
--------array access p.addr (((i * 1024) + j) - 1)
--------array access q.addr ((i * 1024) + j)
------array access u.addr (((i * 1024) + 1024) - 1)
------j
------Loop Bound: (1022, 1)
------Loop inc: (j + -1)
------Loop predicate: >=
--------array access p.addr ((i * 1024) + j)
--------array access u.addr (((i * 1024) + j) + 1)
--------array access q.addr ((i * 1024) + j)
--------array access u.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 1
--------Sample number: 10
--------Sample number: 10
------Sample number: 1
--------Sample number: 10
--------Sample number: 10
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
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
    for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {
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
int calAddrv_addr0( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
int calAddrp_addr0( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddrv_addr1( int t, int i) {
    int result = ((0 + i)) * 8 / 64;
    return result;
}
int calAddrq_addr0( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddrp_addr1( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrp_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddru_addr0( int t, int i, int j) {
    int result = ((((j * 1024) + i) - 1)) * 8 / 64;
    return result;
}
int calAddru_addr1( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
int calAddru_addr2( int t, int i, int j) {
    int result = ((((j * 1024) + i) + 1)) * 8 / 64;
    return result;
}
int calAddrq_addr1( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrp_addr3( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrq_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr2( int t, int i) {
    int result = ((1047552 + i)) * 8 / 64;
    return result;
}
int calAddrp_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr3( int t, int i, int j) {
    int result = ((((j + 1) * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrq_addr3( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr4( int t, int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
int calAddru_addr3( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddrp_addr5( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddru_addr4( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddrq_addr4( int t, int i) {
    int result = (((i * 1024) + 0)) * 8 / 64;
    return result;
}
int calAddrp_addr6( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrp_addr7( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr5( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr6( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrv_addr7( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrq_addr5( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrp_addr8( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
int calAddrq_addr6( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddru_addr5( int t, int i) {
    int result = ((((i * 1024) + 1024) - 1)) * 8 / 64;
    return result;
}
int calAddrp_addr9( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddru_addr6( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
int calAddrq_addr7( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddru_addr7( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_p_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr0 = -1;
    uint64_t prev_t_Start_p_addr0 = -1;
    uint64_t prev_t_End_p_addr0 = -1;
    uint64_t prev_i_Start_p_addr0 = -1;
    uint64_t prev_i_End_p_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr0 != -1) {
            if ( calAddrp_addr0( t_Start - prev_t_Start_p_addr0 + prev_t_End_p_addr0, i_Start - prev_i_Start_p_addr0 + prev_i_End_p_addr0) == calAddrp_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_p_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_p_addr0 = cnt;
                        prev_t_Start_p_addr0 = t_Start;
                        prev_t_End_p_addr0 = t;
                        prev_i_Start_p_addr0 = i_Start;
                        prev_i_End_p_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr1 = -1;
    uint64_t prev_t_Start_p_addr1 = -1;
    uint64_t prev_t_End_p_addr1 = -1;
    uint64_t prev_i_Start_p_addr1 = -1;
    uint64_t prev_i_End_p_addr1 = -1;
    uint64_t prev_j_Start_p_addr1 = -1;
    uint64_t prev_j_End_p_addr1 = -1;
    uint64_t prev_cnt_p_addr2 = -1;
    uint64_t prev_t_Start_p_addr2 = -1;
    uint64_t prev_t_End_p_addr2 = -1;
    uint64_t prev_i_Start_p_addr2 = -1;
    uint64_t prev_i_End_p_addr2 = -1;
    uint64_t prev_j_Start_p_addr2 = -1;
    uint64_t prev_j_End_p_addr2 = -1;
    uint64_t prev_cnt_p_addr3 = -1;
    uint64_t prev_t_Start_p_addr3 = -1;
    uint64_t prev_t_End_p_addr3 = -1;
    uint64_t prev_i_Start_p_addr3 = -1;
    uint64_t prev_i_End_p_addr3 = -1;
    uint64_t prev_j_Start_p_addr3 = -1;
    uint64_t prev_j_End_p_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr1 != -1) {
            if ( calAddrp_addr1( t_Start - prev_t_Start_p_addr1 + prev_t_End_p_addr1, i_Start - prev_i_Start_p_addr1 + prev_i_End_p_addr1, j_Start - prev_j_Start_p_addr1 + prev_j_End_p_addr1) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr2 != -1) {
            if ( calAddrp_addr2( t_Start - prev_t_Start_p_addr2 + prev_t_End_p_addr2, i_Start - prev_i_Start_p_addr2 + prev_i_End_p_addr2, j_Start - prev_j_Start_p_addr2 + prev_j_End_p_addr2) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr3 != -1) {
            if ( calAddrp_addr3( t_Start - prev_t_Start_p_addr3 + prev_t_End_p_addr3, i_Start - prev_i_Start_p_addr3 + prev_i_End_p_addr3, j_Start - prev_j_Start_p_addr3 + prev_j_End_p_addr3) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr1 = cnt;
                            prev_t_Start_p_addr1 = t_Start;
                            prev_t_End_p_addr1 = t;
                            prev_i_Start_p_addr1 = i_Start;
                            prev_i_End_p_addr1 = i;
                            prev_j_Start_p_addr1 = j_Start;
                            prev_j_End_p_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr2 = cnt;
                            prev_t_Start_p_addr2 = t_Start;
                            prev_t_End_p_addr2 = t;
                            prev_i_Start_p_addr2 = i_Start;
                            prev_i_End_p_addr2 = i;
                            prev_j_Start_p_addr2 = j_Start;
                            prev_j_End_p_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr3 = cnt;
                            prev_t_Start_p_addr3 = t_Start;
                            prev_t_End_p_addr3 = t;
                            prev_i_Start_p_addr3 = i_Start;
                            prev_i_End_p_addr3 = i;
                            prev_j_Start_p_addr3 = j_Start;
                            prev_j_End_p_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr1 = -1;
    uint64_t prev_t_Start_p_addr1 = -1;
    uint64_t prev_t_End_p_addr1 = -1;
    uint64_t prev_i_Start_p_addr1 = -1;
    uint64_t prev_i_End_p_addr1 = -1;
    uint64_t prev_j_Start_p_addr1 = -1;
    uint64_t prev_j_End_p_addr1 = -1;
    uint64_t prev_cnt_p_addr2 = -1;
    uint64_t prev_t_Start_p_addr2 = -1;
    uint64_t prev_t_End_p_addr2 = -1;
    uint64_t prev_i_Start_p_addr2 = -1;
    uint64_t prev_i_End_p_addr2 = -1;
    uint64_t prev_j_Start_p_addr2 = -1;
    uint64_t prev_j_End_p_addr2 = -1;
    uint64_t prev_cnt_p_addr3 = -1;
    uint64_t prev_t_Start_p_addr3 = -1;
    uint64_t prev_t_End_p_addr3 = -1;
    uint64_t prev_i_Start_p_addr3 = -1;
    uint64_t prev_i_End_p_addr3 = -1;
    uint64_t prev_j_Start_p_addr3 = -1;
    uint64_t prev_j_End_p_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr1 != -1) {
            if ( calAddrp_addr1( t_Start - prev_t_Start_p_addr1 + prev_t_End_p_addr1, i_Start - prev_i_Start_p_addr1 + prev_i_End_p_addr1, j_Start - prev_j_Start_p_addr1 + prev_j_End_p_addr1) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr2 != -1) {
            if ( calAddrp_addr2( t_Start - prev_t_Start_p_addr2 + prev_t_End_p_addr2, i_Start - prev_i_Start_p_addr2 + prev_i_End_p_addr2, j_Start - prev_j_Start_p_addr2 + prev_j_End_p_addr2) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr3 != -1) {
            if ( calAddrp_addr3( t_Start - prev_t_Start_p_addr3 + prev_t_End_p_addr3, i_Start - prev_i_Start_p_addr3 + prev_i_End_p_addr3, j_Start - prev_j_Start_p_addr3 + prev_j_End_p_addr3) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr1 = cnt;
                            prev_t_Start_p_addr1 = t_Start;
                            prev_t_End_p_addr1 = t;
                            prev_i_Start_p_addr1 = i_Start;
                            prev_i_End_p_addr1 = i;
                            prev_j_Start_p_addr1 = j_Start;
                            prev_j_End_p_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr2 = cnt;
                            prev_t_Start_p_addr2 = t_Start;
                            prev_t_End_p_addr2 = t;
                            prev_i_Start_p_addr2 = i_Start;
                            prev_i_End_p_addr2 = i;
                            prev_j_Start_p_addr2 = j_Start;
                            prev_j_End_p_addr2 = j;
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr3 = cnt;
                            prev_t_Start_p_addr3 = t_Start;
                            prev_t_End_p_addr3 = t;
                            prev_i_Start_p_addr3 = i_Start;
                            prev_i_End_p_addr3 = i;
                            prev_j_Start_p_addr3 = j_Start;
                            prev_j_End_p_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr1 = -1;
    uint64_t prev_t_Start_p_addr1 = -1;
    uint64_t prev_t_End_p_addr1 = -1;
    uint64_t prev_i_Start_p_addr1 = -1;
    uint64_t prev_i_End_p_addr1 = -1;
    uint64_t prev_j_Start_p_addr1 = -1;
    uint64_t prev_j_End_p_addr1 = -1;
    uint64_t prev_cnt_p_addr2 = -1;
    uint64_t prev_t_Start_p_addr2 = -1;
    uint64_t prev_t_End_p_addr2 = -1;
    uint64_t prev_i_Start_p_addr2 = -1;
    uint64_t prev_i_End_p_addr2 = -1;
    uint64_t prev_j_Start_p_addr2 = -1;
    uint64_t prev_j_End_p_addr2 = -1;
    uint64_t prev_cnt_p_addr3 = -1;
    uint64_t prev_t_Start_p_addr3 = -1;
    uint64_t prev_t_End_p_addr3 = -1;
    uint64_t prev_i_Start_p_addr3 = -1;
    uint64_t prev_i_End_p_addr3 = -1;
    uint64_t prev_j_Start_p_addr3 = -1;
    uint64_t prev_j_End_p_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr1 != -1) {
            if ( calAddrp_addr1( t_Start - prev_t_Start_p_addr1 + prev_t_End_p_addr1, i_Start - prev_i_Start_p_addr1 + prev_i_End_p_addr1, j_Start - prev_j_Start_p_addr1 + prev_j_End_p_addr1) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr2 != -1) {
            if ( calAddrp_addr2( t_Start - prev_t_Start_p_addr2 + prev_t_End_p_addr2, i_Start - prev_i_Start_p_addr2 + prev_i_End_p_addr2, j_Start - prev_j_Start_p_addr2 + prev_j_End_p_addr2) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr3 != -1) {
            if ( calAddrp_addr3( t_Start - prev_t_Start_p_addr3 + prev_t_End_p_addr3, i_Start - prev_i_Start_p_addr3 + prev_i_End_p_addr3, j_Start - prev_j_Start_p_addr3 + prev_j_End_p_addr3) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr1 = cnt;
                            prev_t_Start_p_addr1 = t_Start;
                            prev_t_End_p_addr1 = t;
                            prev_i_Start_p_addr1 = i_Start;
                            prev_i_End_p_addr1 = i;
                            prev_j_Start_p_addr1 = j_Start;
                            prev_j_End_p_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr2 = cnt;
                            prev_t_Start_p_addr2 = t_Start;
                            prev_t_End_p_addr2 = t;
                            prev_i_Start_p_addr2 = i_Start;
                            prev_i_End_p_addr2 = i;
                            prev_j_Start_p_addr2 = j_Start;
                            prev_j_End_p_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr3 = cnt;
                            prev_t_Start_p_addr3 = t_Start;
                            prev_t_End_p_addr3 = t;
                            prev_i_Start_p_addr3 = i_Start;
                            prev_i_End_p_addr3 = i;
                            prev_j_Start_p_addr3 = j_Start;
                            prev_j_End_p_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr4 = -1;
    uint64_t prev_t_Start_p_addr4 = -1;
    uint64_t prev_t_End_p_addr4 = -1;
    uint64_t prev_i_Start_p_addr4 = -1;
    uint64_t prev_i_End_p_addr4 = -1;
    uint64_t prev_j_Start_p_addr4 = -1;
    uint64_t prev_j_End_p_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr4 != -1) {
            if ( calAddrp_addr4( t_Start - prev_t_Start_p_addr4 + prev_t_End_p_addr4, i_Start - prev_i_Start_p_addr4 + prev_i_End_p_addr4, j_Start - prev_j_Start_p_addr4 + prev_j_End_p_addr4) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr4 = cnt;
                            prev_t_Start_p_addr4 = t_Start;
                            prev_t_End_p_addr4 = t;
                            prev_i_Start_p_addr4 = i_Start;
                            prev_i_End_p_addr4 = i;
                            prev_j_Start_p_addr4 = j_Start;
                            prev_j_End_p_addr4 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr5 = -1;
    uint64_t prev_t_Start_p_addr5 = -1;
    uint64_t prev_t_End_p_addr5 = -1;
    uint64_t prev_i_Start_p_addr5 = -1;
    uint64_t prev_i_End_p_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr5 != -1) {
            if ( calAddrp_addr5( t_Start - prev_t_Start_p_addr5 + prev_t_End_p_addr5, i_Start - prev_i_Start_p_addr5 + prev_i_End_p_addr5) == calAddrp_addr5(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_p_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr5(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr5(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_p_addr5 = cnt;
                        prev_t_Start_p_addr5 = t_Start;
                        prev_t_End_p_addr5 = t;
                        prev_i_Start_p_addr5 = i_Start;
                        prev_i_End_p_addr5 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr6 = -1;
    uint64_t prev_t_Start_p_addr6 = -1;
    uint64_t prev_t_End_p_addr6 = -1;
    uint64_t prev_i_Start_p_addr6 = -1;
    uint64_t prev_i_End_p_addr6 = -1;
    uint64_t prev_j_Start_p_addr6 = -1;
    uint64_t prev_j_End_p_addr6 = -1;
    uint64_t prev_cnt_p_addr7 = -1;
    uint64_t prev_t_Start_p_addr7 = -1;
    uint64_t prev_t_End_p_addr7 = -1;
    uint64_t prev_i_Start_p_addr7 = -1;
    uint64_t prev_i_End_p_addr7 = -1;
    uint64_t prev_j_Start_p_addr7 = -1;
    uint64_t prev_j_End_p_addr7 = -1;
    uint64_t prev_cnt_p_addr8 = -1;
    uint64_t prev_t_Start_p_addr8 = -1;
    uint64_t prev_t_End_p_addr8 = -1;
    uint64_t prev_i_Start_p_addr8 = -1;
    uint64_t prev_i_End_p_addr8 = -1;
    uint64_t prev_j_Start_p_addr8 = -1;
    uint64_t prev_j_End_p_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr6 != -1) {
            if ( calAddrp_addr6( t_Start - prev_t_Start_p_addr6 + prev_t_End_p_addr6, i_Start - prev_i_Start_p_addr6 + prev_i_End_p_addr6, j_Start - prev_j_Start_p_addr6 + prev_j_End_p_addr6) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr7 != -1) {
            if ( calAddrp_addr7( t_Start - prev_t_Start_p_addr7 + prev_t_End_p_addr7, i_Start - prev_i_Start_p_addr7 + prev_i_End_p_addr7, j_Start - prev_j_Start_p_addr7 + prev_j_End_p_addr7) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr8 != -1) {
            if ( calAddrp_addr8( t_Start - prev_t_Start_p_addr8 + prev_t_End_p_addr8, i_Start - prev_i_Start_p_addr8 + prev_i_End_p_addr8, j_Start - prev_j_Start_p_addr8 + prev_j_End_p_addr8) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr6 = cnt;
                            prev_t_Start_p_addr6 = t_Start;
                            prev_t_End_p_addr6 = t;
                            prev_i_Start_p_addr6 = i_Start;
                            prev_i_End_p_addr6 = i;
                            prev_j_Start_p_addr6 = j_Start;
                            prev_j_End_p_addr6 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr7 = cnt;
                            prev_t_Start_p_addr7 = t_Start;
                            prev_t_End_p_addr7 = t;
                            prev_i_Start_p_addr7 = i_Start;
                            prev_i_End_p_addr7 = i;
                            prev_j_Start_p_addr7 = j_Start;
                            prev_j_End_p_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr8 = cnt;
                            prev_t_Start_p_addr8 = t_Start;
                            prev_t_End_p_addr8 = t;
                            prev_i_Start_p_addr8 = i_Start;
                            prev_i_End_p_addr8 = i;
                            prev_j_Start_p_addr8 = j_Start;
                            prev_j_End_p_addr8 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr6 = -1;
    uint64_t prev_t_Start_p_addr6 = -1;
    uint64_t prev_t_End_p_addr6 = -1;
    uint64_t prev_i_Start_p_addr6 = -1;
    uint64_t prev_i_End_p_addr6 = -1;
    uint64_t prev_j_Start_p_addr6 = -1;
    uint64_t prev_j_End_p_addr6 = -1;
    uint64_t prev_cnt_p_addr7 = -1;
    uint64_t prev_t_Start_p_addr7 = -1;
    uint64_t prev_t_End_p_addr7 = -1;
    uint64_t prev_i_Start_p_addr7 = -1;
    uint64_t prev_i_End_p_addr7 = -1;
    uint64_t prev_j_Start_p_addr7 = -1;
    uint64_t prev_j_End_p_addr7 = -1;
    uint64_t prev_cnt_p_addr8 = -1;
    uint64_t prev_t_Start_p_addr8 = -1;
    uint64_t prev_t_End_p_addr8 = -1;
    uint64_t prev_i_Start_p_addr8 = -1;
    uint64_t prev_i_End_p_addr8 = -1;
    uint64_t prev_j_Start_p_addr8 = -1;
    uint64_t prev_j_End_p_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr6 != -1) {
            if ( calAddrp_addr6( t_Start - prev_t_Start_p_addr6 + prev_t_End_p_addr6, i_Start - prev_i_Start_p_addr6 + prev_i_End_p_addr6, j_Start - prev_j_Start_p_addr6 + prev_j_End_p_addr6) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr7 != -1) {
            if ( calAddrp_addr7( t_Start - prev_t_Start_p_addr7 + prev_t_End_p_addr7, i_Start - prev_i_Start_p_addr7 + prev_i_End_p_addr7, j_Start - prev_j_Start_p_addr7 + prev_j_End_p_addr7) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr8 != -1) {
            if ( calAddrp_addr8( t_Start - prev_t_Start_p_addr8 + prev_t_End_p_addr8, i_Start - prev_i_Start_p_addr8 + prev_i_End_p_addr8, j_Start - prev_j_Start_p_addr8 + prev_j_End_p_addr8) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr6 = cnt;
                            prev_t_Start_p_addr6 = t_Start;
                            prev_t_End_p_addr6 = t;
                            prev_i_Start_p_addr6 = i_Start;
                            prev_i_End_p_addr6 = i;
                            prev_j_Start_p_addr6 = j_Start;
                            prev_j_End_p_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr7 = cnt;
                            prev_t_Start_p_addr7 = t_Start;
                            prev_t_End_p_addr7 = t;
                            prev_i_Start_p_addr7 = i_Start;
                            prev_i_End_p_addr7 = i;
                            prev_j_Start_p_addr7 = j_Start;
                            prev_j_End_p_addr7 = j;
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr8 = cnt;
                            prev_t_Start_p_addr8 = t_Start;
                            prev_t_End_p_addr8 = t;
                            prev_i_Start_p_addr8 = i_Start;
                            prev_i_End_p_addr8 = i;
                            prev_j_Start_p_addr8 = j_Start;
                            prev_j_End_p_addr8 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr6 = -1;
    uint64_t prev_t_Start_p_addr6 = -1;
    uint64_t prev_t_End_p_addr6 = -1;
    uint64_t prev_i_Start_p_addr6 = -1;
    uint64_t prev_i_End_p_addr6 = -1;
    uint64_t prev_j_Start_p_addr6 = -1;
    uint64_t prev_j_End_p_addr6 = -1;
    uint64_t prev_cnt_p_addr7 = -1;
    uint64_t prev_t_Start_p_addr7 = -1;
    uint64_t prev_t_End_p_addr7 = -1;
    uint64_t prev_i_Start_p_addr7 = -1;
    uint64_t prev_i_End_p_addr7 = -1;
    uint64_t prev_j_Start_p_addr7 = -1;
    uint64_t prev_j_End_p_addr7 = -1;
    uint64_t prev_cnt_p_addr8 = -1;
    uint64_t prev_t_Start_p_addr8 = -1;
    uint64_t prev_t_End_p_addr8 = -1;
    uint64_t prev_i_Start_p_addr8 = -1;
    uint64_t prev_i_End_p_addr8 = -1;
    uint64_t prev_j_Start_p_addr8 = -1;
    uint64_t prev_j_End_p_addr8 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr6 != -1) {
            if ( calAddrp_addr6( t_Start - prev_t_Start_p_addr6 + prev_t_End_p_addr6, i_Start - prev_i_Start_p_addr6 + prev_i_End_p_addr6, j_Start - prev_j_Start_p_addr6 + prev_j_End_p_addr6) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr7 != -1) {
            if ( calAddrp_addr7( t_Start - prev_t_Start_p_addr7 + prev_t_End_p_addr7, i_Start - prev_i_Start_p_addr7 + prev_i_End_p_addr7, j_Start - prev_j_Start_p_addr7 + prev_j_End_p_addr7) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_p_addr8 != -1) {
            if ( calAddrp_addr8( t_Start - prev_t_Start_p_addr8 + prev_t_End_p_addr8, i_Start - prev_i_Start_p_addr8 + prev_i_End_p_addr8, j_Start - prev_j_Start_p_addr8 + prev_j_End_p_addr8) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr8);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr6 = cnt;
                            prev_t_Start_p_addr6 = t_Start;
                            prev_t_End_p_addr6 = t;
                            prev_i_Start_p_addr6 = i_Start;
                            prev_i_End_p_addr6 = i;
                            prev_j_Start_p_addr6 = j_Start;
                            prev_j_End_p_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr7 = cnt;
                            prev_t_Start_p_addr7 = t_Start;
                            prev_t_End_p_addr7 = t;
                            prev_i_Start_p_addr7 = i_Start;
                            prev_i_End_p_addr7 = i;
                            prev_j_Start_p_addr7 = j_Start;
                            prev_j_End_p_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr8 = cnt;
                            prev_t_Start_p_addr8 = t_Start;
                            prev_t_End_p_addr8 = t;
                            prev_i_Start_p_addr8 = i_Start;
                            prev_i_End_p_addr8 = i;
                            prev_j_Start_p_addr8 = j_Start;
                            prev_j_End_p_addr8 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr8(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_p_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_p_addr9 = -1;
    uint64_t prev_t_Start_p_addr9 = -1;
    uint64_t prev_t_End_p_addr9 = -1;
    uint64_t prev_i_Start_p_addr9 = -1;
    uint64_t prev_i_End_p_addr9 = -1;
    uint64_t prev_j_Start_p_addr9 = -1;
    uint64_t prev_j_End_p_addr9 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_p_addr9 != -1) {
            if ( calAddrp_addr9( t_Start - prev_t_Start_p_addr9 + prev_t_End_p_addr9, i_Start - prev_i_Start_p_addr9 + prev_i_End_p_addr9, j_Start - prev_j_Start_p_addr9 + prev_j_End_p_addr9) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_p_addr9);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr0( t, i) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr1( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr2( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr3( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr4( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrp_addr5( t, i) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr6( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr7( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
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
                        if ( calAddrp_addr8( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrp_addr9( t, i, j) == calAddrp_addr9(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_p_addr9 = cnt;
                            prev_t_Start_p_addr9 = t_Start;
                            prev_t_End_p_addr9 = t;
                            prev_i_Start_p_addr9 = i_Start;
                            prev_i_End_p_addr9 = i;
                            prev_j_Start_p_addr9 = j_Start;
                            prev_j_End_p_addr9 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr0 = -1;
    uint64_t prev_t_Start_q_addr0 = -1;
    uint64_t prev_t_End_q_addr0 = -1;
    uint64_t prev_i_Start_q_addr0 = -1;
    uint64_t prev_i_End_q_addr0 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr0 != -1) {
            if ( calAddrq_addr0( t_Start - prev_t_Start_q_addr0 + prev_t_End_q_addr0, i_Start - prev_i_Start_q_addr0 + prev_i_End_q_addr0) == calAddrq_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_q_addr0);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr0 = cnt;
                        prev_t_Start_q_addr0 = t_Start;
                        prev_t_End_q_addr0 = t;
                        prev_i_Start_q_addr0 = i_Start;
                        prev_i_End_q_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr1 = -1;
    uint64_t prev_t_Start_q_addr1 = -1;
    uint64_t prev_t_End_q_addr1 = -1;
    uint64_t prev_i_Start_q_addr1 = -1;
    uint64_t prev_i_End_q_addr1 = -1;
    uint64_t prev_j_Start_q_addr1 = -1;
    uint64_t prev_j_End_q_addr1 = -1;
    uint64_t prev_cnt_q_addr2 = -1;
    uint64_t prev_t_Start_q_addr2 = -1;
    uint64_t prev_t_End_q_addr2 = -1;
    uint64_t prev_i_Start_q_addr2 = -1;
    uint64_t prev_i_End_q_addr2 = -1;
    uint64_t prev_j_Start_q_addr2 = -1;
    uint64_t prev_j_End_q_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr1 != -1) {
            if ( calAddrq_addr1( t_Start - prev_t_Start_q_addr1 + prev_t_End_q_addr1, i_Start - prev_i_Start_q_addr1 + prev_i_End_q_addr1, j_Start - prev_j_Start_q_addr1 + prev_j_End_q_addr1) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr2 != -1) {
            if ( calAddrq_addr2( t_Start - prev_t_Start_q_addr2 + prev_t_End_q_addr2, i_Start - prev_i_Start_q_addr2 + prev_i_End_q_addr2, j_Start - prev_j_Start_q_addr2 + prev_j_End_q_addr2) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr1 = cnt;
                            prev_t_Start_q_addr1 = t_Start;
                            prev_t_End_q_addr1 = t;
                            prev_i_Start_q_addr1 = i_Start;
                            prev_i_End_q_addr1 = i;
                            prev_j_Start_q_addr1 = j_Start;
                            prev_j_End_q_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr2 = cnt;
                            prev_t_Start_q_addr2 = t_Start;
                            prev_t_End_q_addr2 = t;
                            prev_i_Start_q_addr2 = i_Start;
                            prev_i_End_q_addr2 = i;
                            prev_j_Start_q_addr2 = j_Start;
                            prev_j_End_q_addr2 = j;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr1 = -1;
    uint64_t prev_t_Start_q_addr1 = -1;
    uint64_t prev_t_End_q_addr1 = -1;
    uint64_t prev_i_Start_q_addr1 = -1;
    uint64_t prev_i_End_q_addr1 = -1;
    uint64_t prev_j_Start_q_addr1 = -1;
    uint64_t prev_j_End_q_addr1 = -1;
    uint64_t prev_cnt_q_addr2 = -1;
    uint64_t prev_t_Start_q_addr2 = -1;
    uint64_t prev_t_End_q_addr2 = -1;
    uint64_t prev_i_Start_q_addr2 = -1;
    uint64_t prev_i_End_q_addr2 = -1;
    uint64_t prev_j_Start_q_addr2 = -1;
    uint64_t prev_j_End_q_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr1 != -1) {
            if ( calAddrq_addr1( t_Start - prev_t_Start_q_addr1 + prev_t_End_q_addr1, i_Start - prev_i_Start_q_addr1 + prev_i_End_q_addr1, j_Start - prev_j_Start_q_addr1 + prev_j_End_q_addr1) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr2 != -1) {
            if ( calAddrq_addr2( t_Start - prev_t_Start_q_addr2 + prev_t_End_q_addr2, i_Start - prev_i_Start_q_addr2 + prev_i_End_q_addr2, j_Start - prev_j_Start_q_addr2 + prev_j_End_q_addr2) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr1 = cnt;
                            prev_t_Start_q_addr1 = t_Start;
                            prev_t_End_q_addr1 = t;
                            prev_i_Start_q_addr1 = i_Start;
                            prev_i_End_q_addr1 = i;
                            prev_j_Start_q_addr1 = j_Start;
                            prev_j_End_q_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr2 = cnt;
                            prev_t_Start_q_addr2 = t_Start;
                            prev_t_End_q_addr2 = t;
                            prev_i_Start_q_addr2 = i_Start;
                            prev_i_End_q_addr2 = i;
                            prev_j_Start_q_addr2 = j_Start;
                            prev_j_End_q_addr2 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr3 = -1;
    uint64_t prev_t_Start_q_addr3 = -1;
    uint64_t prev_t_End_q_addr3 = -1;
    uint64_t prev_i_Start_q_addr3 = -1;
    uint64_t prev_i_End_q_addr3 = -1;
    uint64_t prev_j_Start_q_addr3 = -1;
    uint64_t prev_j_End_q_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr3 != -1) {
            if ( calAddrq_addr3( t_Start - prev_t_Start_q_addr3 + prev_t_End_q_addr3, i_Start - prev_i_Start_q_addr3 + prev_i_End_q_addr3, j_Start - prev_j_Start_q_addr3 + prev_j_End_q_addr3) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr3 = cnt;
                            prev_t_Start_q_addr3 = t_Start;
                            prev_t_End_q_addr3 = t;
                            prev_i_Start_q_addr3 = i_Start;
                            prev_i_End_q_addr3 = i;
                            prev_j_Start_q_addr3 = j_Start;
                            prev_j_End_q_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr4 = -1;
    uint64_t prev_t_Start_q_addr4 = -1;
    uint64_t prev_t_End_q_addr4 = -1;
    uint64_t prev_i_Start_q_addr4 = -1;
    uint64_t prev_i_End_q_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr4 != -1) {
            if ( calAddrq_addr4( t_Start - prev_t_Start_q_addr4 + prev_t_End_q_addr4, i_Start - prev_i_Start_q_addr4 + prev_i_End_q_addr4) == calAddrq_addr4(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_q_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr4(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr4(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_q_addr4 = cnt;
                        prev_t_Start_q_addr4 = t_Start;
                        prev_t_End_q_addr4 = t;
                        prev_i_Start_q_addr4 = i_Start;
                        prev_i_End_q_addr4 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr5 = -1;
    uint64_t prev_t_Start_q_addr5 = -1;
    uint64_t prev_t_End_q_addr5 = -1;
    uint64_t prev_i_Start_q_addr5 = -1;
    uint64_t prev_i_End_q_addr5 = -1;
    uint64_t prev_j_Start_q_addr5 = -1;
    uint64_t prev_j_End_q_addr5 = -1;
    uint64_t prev_cnt_q_addr6 = -1;
    uint64_t prev_t_Start_q_addr6 = -1;
    uint64_t prev_t_End_q_addr6 = -1;
    uint64_t prev_i_Start_q_addr6 = -1;
    uint64_t prev_i_End_q_addr6 = -1;
    uint64_t prev_j_Start_q_addr6 = -1;
    uint64_t prev_j_End_q_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr5 != -1) {
            if ( calAddrq_addr5( t_Start - prev_t_Start_q_addr5 + prev_t_End_q_addr5, i_Start - prev_i_Start_q_addr5 + prev_i_End_q_addr5, j_Start - prev_j_Start_q_addr5 + prev_j_End_q_addr5) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr6 != -1) {
            if ( calAddrq_addr6( t_Start - prev_t_Start_q_addr6 + prev_t_End_q_addr6, i_Start - prev_i_Start_q_addr6 + prev_i_End_q_addr6, j_Start - prev_j_Start_q_addr6 + prev_j_End_q_addr6) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr5 = cnt;
                            prev_t_Start_q_addr5 = t_Start;
                            prev_t_End_q_addr5 = t;
                            prev_i_Start_q_addr5 = i_Start;
                            prev_i_End_q_addr5 = i;
                            prev_j_Start_q_addr5 = j_Start;
                            prev_j_End_q_addr5 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr6 = cnt;
                            prev_t_Start_q_addr6 = t_Start;
                            prev_t_End_q_addr6 = t;
                            prev_i_Start_q_addr6 = i_Start;
                            prev_i_End_q_addr6 = i;
                            prev_j_Start_q_addr6 = j_Start;
                            prev_j_End_q_addr6 = j;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr5 = -1;
    uint64_t prev_t_Start_q_addr5 = -1;
    uint64_t prev_t_End_q_addr5 = -1;
    uint64_t prev_i_Start_q_addr5 = -1;
    uint64_t prev_i_End_q_addr5 = -1;
    uint64_t prev_j_Start_q_addr5 = -1;
    uint64_t prev_j_End_q_addr5 = -1;
    uint64_t prev_cnt_q_addr6 = -1;
    uint64_t prev_t_Start_q_addr6 = -1;
    uint64_t prev_t_End_q_addr6 = -1;
    uint64_t prev_i_Start_q_addr6 = -1;
    uint64_t prev_i_End_q_addr6 = -1;
    uint64_t prev_j_Start_q_addr6 = -1;
    uint64_t prev_j_End_q_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr5 != -1) {
            if ( calAddrq_addr5( t_Start - prev_t_Start_q_addr5 + prev_t_End_q_addr5, i_Start - prev_i_Start_q_addr5 + prev_i_End_q_addr5, j_Start - prev_j_Start_q_addr5 + prev_j_End_q_addr5) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_q_addr6 != -1) {
            if ( calAddrq_addr6( t_Start - prev_t_Start_q_addr6 + prev_t_End_q_addr6, i_Start - prev_i_Start_q_addr6 + prev_i_End_q_addr6, j_Start - prev_j_Start_q_addr6 + prev_j_End_q_addr6) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr5 = cnt;
                            prev_t_Start_q_addr5 = t_Start;
                            prev_t_End_q_addr5 = t;
                            prev_i_Start_q_addr5 = i_Start;
                            prev_i_End_q_addr5 = i;
                            prev_j_Start_q_addr5 = j_Start;
                            prev_j_End_q_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr6 = cnt;
                            prev_t_Start_q_addr6 = t_Start;
                            prev_t_End_q_addr6 = t;
                            prev_i_Start_q_addr6 = i_Start;
                            prev_i_End_q_addr6 = i;
                            prev_j_Start_q_addr6 = j_Start;
                            prev_j_End_q_addr6 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_q_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_q_addr7 = -1;
    uint64_t prev_t_Start_q_addr7 = -1;
    uint64_t prev_t_End_q_addr7 = -1;
    uint64_t prev_i_Start_q_addr7 = -1;
    uint64_t prev_i_End_q_addr7 = -1;
    uint64_t prev_j_Start_q_addr7 = -1;
    uint64_t prev_j_End_q_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_q_addr7 != -1) {
            if ( calAddrq_addr7( t_Start - prev_t_Start_q_addr7 + prev_t_End_q_addr7, i_Start - prev_i_Start_q_addr7 + prev_i_End_q_addr7, j_Start - prev_j_Start_q_addr7 + prev_j_End_q_addr7) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_q_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr0( t, i) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr1( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr2( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr3( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrq_addr4( t, i) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr5( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr6( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrq_addr7( t, i, j) == calAddrq_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_q_addr7 = cnt;
                            prev_t_Start_q_addr7 = t_Start;
                            prev_t_End_q_addr7 = t;
                            prev_i_Start_q_addr7 = i_Start;
                            prev_i_End_q_addr7 = i;
                            prev_j_Start_q_addr7 = j_Start;
                            prev_j_End_q_addr7 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr0 = -1;
    uint64_t prev_t_Start_u_addr0 = -1;
    uint64_t prev_t_End_u_addr0 = -1;
    uint64_t prev_i_Start_u_addr0 = -1;
    uint64_t prev_i_End_u_addr0 = -1;
    uint64_t prev_j_Start_u_addr0 = -1;
    uint64_t prev_j_End_u_addr0 = -1;
    uint64_t prev_cnt_u_addr1 = -1;
    uint64_t prev_t_Start_u_addr1 = -1;
    uint64_t prev_t_End_u_addr1 = -1;
    uint64_t prev_i_Start_u_addr1 = -1;
    uint64_t prev_i_End_u_addr1 = -1;
    uint64_t prev_j_Start_u_addr1 = -1;
    uint64_t prev_j_End_u_addr1 = -1;
    uint64_t prev_cnt_u_addr2 = -1;
    uint64_t prev_t_Start_u_addr2 = -1;
    uint64_t prev_t_End_u_addr2 = -1;
    uint64_t prev_i_Start_u_addr2 = -1;
    uint64_t prev_i_End_u_addr2 = -1;
    uint64_t prev_j_Start_u_addr2 = -1;
    uint64_t prev_j_End_u_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr0 != -1) {
            if ( calAddru_addr0( t_Start - prev_t_Start_u_addr0 + prev_t_End_u_addr0, i_Start - prev_i_Start_u_addr0 + prev_i_End_u_addr0, j_Start - prev_j_Start_u_addr0 + prev_j_End_u_addr0) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr1 != -1) {
            if ( calAddru_addr1( t_Start - prev_t_Start_u_addr1 + prev_t_End_u_addr1, i_Start - prev_i_Start_u_addr1 + prev_i_End_u_addr1, j_Start - prev_j_Start_u_addr1 + prev_j_End_u_addr1) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr2 != -1) {
            if ( calAddru_addr2( t_Start - prev_t_Start_u_addr2 + prev_t_End_u_addr2, i_Start - prev_i_Start_u_addr2 + prev_i_End_u_addr2, j_Start - prev_j_Start_u_addr2 + prev_j_End_u_addr2) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr0 = cnt;
                            prev_t_Start_u_addr0 = t_Start;
                            prev_t_End_u_addr0 = t;
                            prev_i_Start_u_addr0 = i_Start;
                            prev_i_End_u_addr0 = i;
                            prev_j_Start_u_addr0 = j_Start;
                            prev_j_End_u_addr0 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr1 = cnt;
                            prev_t_Start_u_addr1 = t_Start;
                            prev_t_End_u_addr1 = t;
                            prev_i_Start_u_addr1 = i_Start;
                            prev_i_End_u_addr1 = i;
                            prev_j_Start_u_addr1 = j_Start;
                            prev_j_End_u_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr2 = cnt;
                            prev_t_Start_u_addr2 = t_Start;
                            prev_t_End_u_addr2 = t;
                            prev_i_Start_u_addr2 = i_Start;
                            prev_i_End_u_addr2 = i;
                            prev_j_Start_u_addr2 = j_Start;
                            prev_j_End_u_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr0(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr0 = -1;
    uint64_t prev_t_Start_u_addr0 = -1;
    uint64_t prev_t_End_u_addr0 = -1;
    uint64_t prev_i_Start_u_addr0 = -1;
    uint64_t prev_i_End_u_addr0 = -1;
    uint64_t prev_j_Start_u_addr0 = -1;
    uint64_t prev_j_End_u_addr0 = -1;
    uint64_t prev_cnt_u_addr1 = -1;
    uint64_t prev_t_Start_u_addr1 = -1;
    uint64_t prev_t_End_u_addr1 = -1;
    uint64_t prev_i_Start_u_addr1 = -1;
    uint64_t prev_i_End_u_addr1 = -1;
    uint64_t prev_j_Start_u_addr1 = -1;
    uint64_t prev_j_End_u_addr1 = -1;
    uint64_t prev_cnt_u_addr2 = -1;
    uint64_t prev_t_Start_u_addr2 = -1;
    uint64_t prev_t_End_u_addr2 = -1;
    uint64_t prev_i_Start_u_addr2 = -1;
    uint64_t prev_i_End_u_addr2 = -1;
    uint64_t prev_j_Start_u_addr2 = -1;
    uint64_t prev_j_End_u_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr0 != -1) {
            if ( calAddru_addr0( t_Start - prev_t_Start_u_addr0 + prev_t_End_u_addr0, i_Start - prev_i_Start_u_addr0 + prev_i_End_u_addr0, j_Start - prev_j_Start_u_addr0 + prev_j_End_u_addr0) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr1 != -1) {
            if ( calAddru_addr1( t_Start - prev_t_Start_u_addr1 + prev_t_End_u_addr1, i_Start - prev_i_Start_u_addr1 + prev_i_End_u_addr1, j_Start - prev_j_Start_u_addr1 + prev_j_End_u_addr1) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr2 != -1) {
            if ( calAddru_addr2( t_Start - prev_t_Start_u_addr2 + prev_t_End_u_addr2, i_Start - prev_i_Start_u_addr2 + prev_i_End_u_addr2, j_Start - prev_j_Start_u_addr2 + prev_j_End_u_addr2) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr0 = cnt;
                            prev_t_Start_u_addr0 = t_Start;
                            prev_t_End_u_addr0 = t;
                            prev_i_Start_u_addr0 = i_Start;
                            prev_i_End_u_addr0 = i;
                            prev_j_Start_u_addr0 = j_Start;
                            prev_j_End_u_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr1 = cnt;
                            prev_t_Start_u_addr1 = t_Start;
                            prev_t_End_u_addr1 = t;
                            prev_i_Start_u_addr1 = i_Start;
                            prev_i_End_u_addr1 = i;
                            prev_j_Start_u_addr1 = j_Start;
                            prev_j_End_u_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr2 = cnt;
                            prev_t_Start_u_addr2 = t_Start;
                            prev_t_End_u_addr2 = t;
                            prev_i_Start_u_addr2 = i_Start;
                            prev_i_End_u_addr2 = i;
                            prev_j_Start_u_addr2 = j_Start;
                            prev_j_End_u_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr1(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr0 = -1;
    uint64_t prev_t_Start_u_addr0 = -1;
    uint64_t prev_t_End_u_addr0 = -1;
    uint64_t prev_i_Start_u_addr0 = -1;
    uint64_t prev_i_End_u_addr0 = -1;
    uint64_t prev_j_Start_u_addr0 = -1;
    uint64_t prev_j_End_u_addr0 = -1;
    uint64_t prev_cnt_u_addr1 = -1;
    uint64_t prev_t_Start_u_addr1 = -1;
    uint64_t prev_t_End_u_addr1 = -1;
    uint64_t prev_i_Start_u_addr1 = -1;
    uint64_t prev_i_End_u_addr1 = -1;
    uint64_t prev_j_Start_u_addr1 = -1;
    uint64_t prev_j_End_u_addr1 = -1;
    uint64_t prev_cnt_u_addr2 = -1;
    uint64_t prev_t_Start_u_addr2 = -1;
    uint64_t prev_t_End_u_addr2 = -1;
    uint64_t prev_i_Start_u_addr2 = -1;
    uint64_t prev_i_End_u_addr2 = -1;
    uint64_t prev_j_Start_u_addr2 = -1;
    uint64_t prev_j_End_u_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr0 != -1) {
            if ( calAddru_addr0( t_Start - prev_t_Start_u_addr0 + prev_t_End_u_addr0, i_Start - prev_i_Start_u_addr0 + prev_i_End_u_addr0, j_Start - prev_j_Start_u_addr0 + prev_j_End_u_addr0) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr1 != -1) {
            if ( calAddru_addr1( t_Start - prev_t_Start_u_addr1 + prev_t_End_u_addr1, i_Start - prev_i_Start_u_addr1 + prev_i_End_u_addr1, j_Start - prev_j_Start_u_addr1 + prev_j_End_u_addr1) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr2 != -1) {
            if ( calAddru_addr2( t_Start - prev_t_Start_u_addr2 + prev_t_End_u_addr2, i_Start - prev_i_Start_u_addr2 + prev_i_End_u_addr2, j_Start - prev_j_Start_u_addr2 + prev_j_End_u_addr2) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr0 = cnt;
                            prev_t_Start_u_addr0 = t_Start;
                            prev_t_End_u_addr0 = t;
                            prev_i_Start_u_addr0 = i_Start;
                            prev_i_End_u_addr0 = i;
                            prev_j_Start_u_addr0 = j_Start;
                            prev_j_End_u_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr1 = cnt;
                            prev_t_Start_u_addr1 = t_Start;
                            prev_t_End_u_addr1 = t;
                            prev_i_Start_u_addr1 = i_Start;
                            prev_i_End_u_addr1 = i;
                            prev_j_Start_u_addr1 = j_Start;
                            prev_j_End_u_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr2 = cnt;
                            prev_t_Start_u_addr2 = t_Start;
                            prev_t_End_u_addr2 = t;
                            prev_i_Start_u_addr2 = i_Start;
                            prev_i_End_u_addr2 = i;
                            prev_j_Start_u_addr2 = j_Start;
                            prev_j_End_u_addr2 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr2(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr3 = -1;
    uint64_t prev_t_Start_u_addr3 = -1;
    uint64_t prev_t_End_u_addr3 = -1;
    uint64_t prev_i_Start_u_addr3 = -1;
    uint64_t prev_i_End_u_addr3 = -1;
    uint64_t prev_cnt_u_addr4 = -1;
    uint64_t prev_t_Start_u_addr4 = -1;
    uint64_t prev_t_End_u_addr4 = -1;
    uint64_t prev_i_Start_u_addr4 = -1;
    uint64_t prev_i_End_u_addr4 = -1;
    uint64_t prev_cnt_u_addr5 = -1;
    uint64_t prev_t_Start_u_addr5 = -1;
    uint64_t prev_t_End_u_addr5 = -1;
    uint64_t prev_i_Start_u_addr5 = -1;
    uint64_t prev_i_End_u_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr3 != -1) {
            if ( calAddru_addr3( t_Start - prev_t_Start_u_addr3 + prev_t_End_u_addr3, i_Start - prev_i_Start_u_addr3 + prev_i_End_u_addr3) == calAddru_addr3(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr4 != -1) {
            if ( calAddru_addr4( t_Start - prev_t_Start_u_addr4 + prev_t_End_u_addr4, i_Start - prev_i_Start_u_addr4 + prev_i_End_u_addr4) == calAddru_addr3(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr5 != -1) {
            if ( calAddru_addr5( t_Start - prev_t_Start_u_addr5 + prev_t_End_u_addr5, i_Start - prev_i_Start_u_addr5 + prev_i_End_u_addr5) == calAddru_addr3(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr3 = cnt;
                        prev_t_Start_u_addr3 = t_Start;
                        prev_t_End_u_addr3 = t;
                        prev_i_Start_u_addr3 = i_Start;
                        prev_i_End_u_addr3 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr4 = cnt;
                        prev_t_Start_u_addr4 = t_Start;
                        prev_t_End_u_addr4 = t;
                        prev_i_Start_u_addr4 = i_Start;
                        prev_i_End_u_addr4 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr3(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr5 = cnt;
                        prev_t_Start_u_addr5 = t_Start;
                        prev_t_End_u_addr5 = t;
                        prev_i_Start_u_addr5 = i_Start;
                        prev_i_End_u_addr5 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr3(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr3 = -1;
    uint64_t prev_t_Start_u_addr3 = -1;
    uint64_t prev_t_End_u_addr3 = -1;
    uint64_t prev_i_Start_u_addr3 = -1;
    uint64_t prev_i_End_u_addr3 = -1;
    uint64_t prev_cnt_u_addr4 = -1;
    uint64_t prev_t_Start_u_addr4 = -1;
    uint64_t prev_t_End_u_addr4 = -1;
    uint64_t prev_i_Start_u_addr4 = -1;
    uint64_t prev_i_End_u_addr4 = -1;
    uint64_t prev_cnt_u_addr5 = -1;
    uint64_t prev_t_Start_u_addr5 = -1;
    uint64_t prev_t_End_u_addr5 = -1;
    uint64_t prev_i_Start_u_addr5 = -1;
    uint64_t prev_i_End_u_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr3 != -1) {
            if ( calAddru_addr3( t_Start - prev_t_Start_u_addr3 + prev_t_End_u_addr3, i_Start - prev_i_Start_u_addr3 + prev_i_End_u_addr3) == calAddru_addr4(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr4 != -1) {
            if ( calAddru_addr4( t_Start - prev_t_Start_u_addr4 + prev_t_End_u_addr4, i_Start - prev_i_Start_u_addr4 + prev_i_End_u_addr4) == calAddru_addr4(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr5 != -1) {
            if ( calAddru_addr5( t_Start - prev_t_Start_u_addr5 + prev_t_End_u_addr5, i_Start - prev_i_Start_u_addr5 + prev_i_End_u_addr5) == calAddru_addr4(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr4(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr3 = cnt;
                        prev_t_Start_u_addr3 = t_Start;
                        prev_t_End_u_addr3 = t;
                        prev_i_Start_u_addr3 = i_Start;
                        prev_i_End_u_addr3 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr4(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr4 = cnt;
                        prev_t_Start_u_addr4 = t_Start;
                        prev_t_End_u_addr4 = t;
                        prev_i_Start_u_addr4 = i_Start;
                        prev_i_End_u_addr4 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr4(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr5 = cnt;
                        prev_t_Start_u_addr5 = t_Start;
                        prev_t_End_u_addr5 = t;
                        prev_i_Start_u_addr5 = i_Start;
                        prev_i_End_u_addr5 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr4(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr3 = -1;
    uint64_t prev_t_Start_u_addr3 = -1;
    uint64_t prev_t_End_u_addr3 = -1;
    uint64_t prev_i_Start_u_addr3 = -1;
    uint64_t prev_i_End_u_addr3 = -1;
    uint64_t prev_cnt_u_addr4 = -1;
    uint64_t prev_t_Start_u_addr4 = -1;
    uint64_t prev_t_End_u_addr4 = -1;
    uint64_t prev_i_Start_u_addr4 = -1;
    uint64_t prev_i_End_u_addr4 = -1;
    uint64_t prev_cnt_u_addr5 = -1;
    uint64_t prev_t_Start_u_addr5 = -1;
    uint64_t prev_t_End_u_addr5 = -1;
    uint64_t prev_i_Start_u_addr5 = -1;
    uint64_t prev_i_End_u_addr5 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr3 != -1) {
            if ( calAddru_addr3( t_Start - prev_t_Start_u_addr3 + prev_t_End_u_addr3, i_Start - prev_i_Start_u_addr3 + prev_i_End_u_addr3) == calAddru_addr5(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr4 != -1) {
            if ( calAddru_addr4( t_Start - prev_t_Start_u_addr4 + prev_t_End_u_addr4, i_Start - prev_i_Start_u_addr4 + prev_i_End_u_addr4) == calAddru_addr5(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr5 != -1) {
            if ( calAddru_addr5( t_Start - prev_t_Start_u_addr5 + prev_t_End_u_addr5, i_Start - prev_i_Start_u_addr5 + prev_i_End_u_addr5) == calAddru_addr5(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_u_addr5);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr5(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr3 = cnt;
                        prev_t_Start_u_addr3 = t_Start;
                        prev_t_End_u_addr3 = t;
                        prev_i_Start_u_addr3 = i_Start;
                        prev_i_End_u_addr3 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr5(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr4 = cnt;
                        prev_t_Start_u_addr4 = t_Start;
                        prev_t_End_u_addr4 = t;
                        prev_i_Start_u_addr4 = i_Start;
                        prev_i_End_u_addr4 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr5(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_u_addr5 = cnt;
                        prev_t_Start_u_addr5 = t_Start;
                        prev_t_End_u_addr5 = t;
                        prev_i_Start_u_addr5 = i_Start;
                        prev_i_End_u_addr5 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr5(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr6 = -1;
    uint64_t prev_t_Start_u_addr6 = -1;
    uint64_t prev_t_End_u_addr6 = -1;
    uint64_t prev_i_Start_u_addr6 = -1;
    uint64_t prev_i_End_u_addr6 = -1;
    uint64_t prev_j_Start_u_addr6 = -1;
    uint64_t prev_j_End_u_addr6 = -1;
    uint64_t prev_cnt_u_addr7 = -1;
    uint64_t prev_t_Start_u_addr7 = -1;
    uint64_t prev_t_End_u_addr7 = -1;
    uint64_t prev_i_Start_u_addr7 = -1;
    uint64_t prev_i_End_u_addr7 = -1;
    uint64_t prev_j_Start_u_addr7 = -1;
    uint64_t prev_j_End_u_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr6 != -1) {
            if ( calAddru_addr6( t_Start - prev_t_Start_u_addr6 + prev_t_End_u_addr6, i_Start - prev_i_Start_u_addr6 + prev_i_End_u_addr6, j_Start - prev_j_Start_u_addr6 + prev_j_End_u_addr6) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr7 != -1) {
            if ( calAddru_addr7( t_Start - prev_t_Start_u_addr7 + prev_t_End_u_addr7, i_Start - prev_i_Start_u_addr7 + prev_i_End_u_addr7, j_Start - prev_j_Start_u_addr7 + prev_j_End_u_addr7) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr6 = cnt;
                            prev_t_Start_u_addr6 = t_Start;
                            prev_t_End_u_addr6 = t;
                            prev_i_Start_u_addr6 = i_Start;
                            prev_i_End_u_addr6 = i;
                            prev_j_Start_u_addr6 = j_Start;
                            prev_j_End_u_addr6 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr7 = cnt;
                            prev_t_Start_u_addr7 = t_Start;
                            prev_t_End_u_addr7 = t;
                            prev_i_Start_u_addr7 = i_Start;
                            prev_i_End_u_addr7 = i;
                            prev_j_Start_u_addr7 = j_Start;
                            prev_j_End_u_addr7 = j;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_u_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_u_addr6 = -1;
    uint64_t prev_t_Start_u_addr6 = -1;
    uint64_t prev_t_End_u_addr6 = -1;
    uint64_t prev_i_Start_u_addr6 = -1;
    uint64_t prev_i_End_u_addr6 = -1;
    uint64_t prev_j_Start_u_addr6 = -1;
    uint64_t prev_j_End_u_addr6 = -1;
    uint64_t prev_cnt_u_addr7 = -1;
    uint64_t prev_t_Start_u_addr7 = -1;
    uint64_t prev_t_End_u_addr7 = -1;
    uint64_t prev_i_Start_u_addr7 = -1;
    uint64_t prev_i_End_u_addr7 = -1;
    uint64_t prev_j_Start_u_addr7 = -1;
    uint64_t prev_j_End_u_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_u_addr6 != -1) {
            if ( calAddru_addr6( t_Start - prev_t_Start_u_addr6 + prev_t_End_u_addr6, i_Start - prev_i_Start_u_addr6 + prev_i_End_u_addr6, j_Start - prev_j_Start_u_addr6 + prev_j_End_u_addr6) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_u_addr7 != -1) {
            if ( calAddru_addr7( t_Start - prev_t_Start_u_addr7 + prev_t_End_u_addr7, i_Start - prev_i_Start_u_addr7 + prev_i_End_u_addr7, j_Start - prev_j_Start_u_addr7 + prev_j_End_u_addr7) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_u_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr0( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr1( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr2( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr3( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr4( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru_addr5( t, i) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB6 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB6 = j_Start;
                }
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr6( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr6 = cnt;
                            prev_t_Start_u_addr6 = t_Start;
                            prev_t_End_u_addr6 = t;
                            prev_i_Start_u_addr6 = i_Start;
                            prev_i_End_u_addr6 = i;
                            prev_j_Start_u_addr6 = j_Start;
                            prev_j_End_u_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddru_addr7( t, i, j) == calAddru_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_u_addr7 = cnt;
                            prev_t_Start_u_addr7 = t_Start;
                            prev_t_End_u_addr7 = t;
                            prev_i_Start_u_addr7 = i_Start;
                            prev_i_End_u_addr7 = i;
                            prev_j_Start_u_addr7 = j_Start;
                            prev_j_End_u_addr7 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr1 = -1;
    uint64_t prev_t_Start_v_addr1 = -1;
    uint64_t prev_t_End_v_addr1 = -1;
    uint64_t prev_i_Start_v_addr1 = -1;
    uint64_t prev_i_End_v_addr1 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr1 != -1) {
            if ( calAddrv_addr1( t_Start - prev_t_Start_v_addr1 + prev_t_End_v_addr1, i_Start - prev_i_Start_v_addr1 + prev_i_End_v_addr1) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr0(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr1 = cnt;
                        prev_t_Start_v_addr1 = t_Start;
                        prev_t_End_v_addr1 = t;
                        prev_i_Start_v_addr1 = i_Start;
                        prev_i_End_v_addr1 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr0(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr0(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr1 = -1;
    uint64_t prev_t_Start_v_addr1 = -1;
    uint64_t prev_t_End_v_addr1 = -1;
    uint64_t prev_i_Start_v_addr1 = -1;
    uint64_t prev_i_End_v_addr1 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr1(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr1 != -1) {
            if ( calAddrv_addr1( t_Start - prev_t_Start_v_addr1 + prev_t_End_v_addr1, i_Start - prev_i_Start_v_addr1 + prev_i_End_v_addr1) == calAddrv_addr1(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr1(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr1 = cnt;
                        prev_t_Start_v_addr1 = t_Start;
                        prev_t_End_v_addr1 = t;
                        prev_i_Start_v_addr1 = i_Start;
                        prev_i_End_v_addr1 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr1(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr1(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr0 = -1;
    uint64_t prev_t_Start_v_addr0 = -1;
    uint64_t prev_t_End_v_addr0 = -1;
    uint64_t prev_i_Start_v_addr0 = -1;
    uint64_t prev_i_End_v_addr0 = -1;
    uint64_t prev_cnt_v_addr1 = -1;
    uint64_t prev_t_Start_v_addr1 = -1;
    uint64_t prev_t_End_v_addr1 = -1;
    uint64_t prev_i_Start_v_addr1 = -1;
    uint64_t prev_i_End_v_addr1 = -1;
    uint64_t prev_cnt_v_addr2 = -1;
    uint64_t prev_t_Start_v_addr2 = -1;
    uint64_t prev_t_End_v_addr2 = -1;
    uint64_t prev_i_Start_v_addr2 = -1;
    uint64_t prev_i_End_v_addr2 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr0 != -1) {
            if ( calAddrv_addr0( t_Start - prev_t_Start_v_addr0 + prev_t_End_v_addr0, i_Start - prev_i_Start_v_addr0 + prev_i_End_v_addr0) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr1 != -1) {
            if ( calAddrv_addr1( t_Start - prev_t_Start_v_addr1 + prev_t_End_v_addr1, i_Start - prev_i_Start_v_addr1 + prev_i_End_v_addr1) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr2 != -1) {
            if ( calAddrv_addr2( t_Start - prev_t_Start_v_addr2 + prev_t_End_v_addr2, i_Start - prev_i_Start_v_addr2 + prev_i_End_v_addr2) == calAddrv_addr2(t_Start, i_Start)) {
                rtHistoCal(prev_cnt_v_addr2);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr0 = cnt;
                        prev_t_Start_v_addr0 = t_Start;
                        prev_t_End_v_addr0 = t;
                        prev_i_Start_v_addr0 = i_Start;
                        prev_i_End_v_addr0 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr1 = cnt;
                        prev_t_Start_v_addr1 = t_Start;
                        prev_t_End_v_addr1 = t;
                        prev_i_Start_v_addr1 = i_Start;
                        prev_i_End_v_addr1 = i;
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr2(t_Start, i_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_v_addr2 = cnt;
                        prev_t_Start_v_addr2 = t_Start;
                        prev_t_End_v_addr2 = t;
                        prev_i_Start_v_addr2 = i_Start;
                        prev_i_End_v_addr2 = i;
                        goto EndSample;
                    }
                }
                cntStart = true;
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr2(t_Start, i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr3 = -1;
    uint64_t prev_t_Start_v_addr3 = -1;
    uint64_t prev_t_End_v_addr3 = -1;
    uint64_t prev_i_Start_v_addr3 = -1;
    uint64_t prev_i_End_v_addr3 = -1;
    uint64_t prev_j_Start_v_addr3 = -1;
    uint64_t prev_j_End_v_addr3 = -1;
    uint64_t prev_cnt_v_addr4 = -1;
    uint64_t prev_t_Start_v_addr4 = -1;
    uint64_t prev_t_End_v_addr4 = -1;
    uint64_t prev_i_Start_v_addr4 = -1;
    uint64_t prev_i_End_v_addr4 = -1;
    uint64_t prev_j_Start_v_addr4 = -1;
    uint64_t prev_j_End_v_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr3 != -1) {
            if ( calAddrv_addr3( t_Start - prev_t_Start_v_addr3 + prev_t_End_v_addr3, i_Start - prev_i_Start_v_addr3 + prev_i_End_v_addr3, j_Start - prev_j_Start_v_addr3 + prev_j_End_v_addr3) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr4 != -1) {
            if ( calAddrv_addr4( t_Start - prev_t_Start_v_addr4 + prev_t_End_v_addr4, i_Start - prev_i_Start_v_addr4 + prev_i_End_v_addr4, j_Start - prev_j_Start_v_addr4 + prev_j_End_v_addr4) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr3 = cnt;
                            prev_t_Start_v_addr3 = t_Start;
                            prev_t_End_v_addr3 = t;
                            prev_i_Start_v_addr3 = i_Start;
                            prev_i_End_v_addr3 = i;
                            prev_j_Start_v_addr3 = j_Start;
                            prev_j_End_v_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr4 = cnt;
                            prev_t_Start_v_addr4 = t_Start;
                            prev_t_End_v_addr4 = t;
                            prev_i_Start_v_addr4 = i_Start;
                            prev_i_End_v_addr4 = i;
                            prev_j_Start_v_addr4 = j_Start;
                            prev_j_End_v_addr4 = j;
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr3(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr3 = -1;
    uint64_t prev_t_Start_v_addr3 = -1;
    uint64_t prev_t_End_v_addr3 = -1;
    uint64_t prev_i_Start_v_addr3 = -1;
    uint64_t prev_i_End_v_addr3 = -1;
    uint64_t prev_j_Start_v_addr3 = -1;
    uint64_t prev_j_End_v_addr3 = -1;
    uint64_t prev_cnt_v_addr4 = -1;
    uint64_t prev_t_Start_v_addr4 = -1;
    uint64_t prev_t_End_v_addr4 = -1;
    uint64_t prev_i_Start_v_addr4 = -1;
    uint64_t prev_i_End_v_addr4 = -1;
    uint64_t prev_j_Start_v_addr4 = -1;
    uint64_t prev_j_End_v_addr4 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1 - 1022 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1 - 1022 + 1) + 1022;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr3 != -1) {
            if ( calAddrv_addr3( t_Start - prev_t_Start_v_addr3 + prev_t_End_v_addr3, i_Start - prev_i_Start_v_addr3 + prev_i_End_v_addr3, j_Start - prev_j_Start_v_addr3 + prev_j_End_v_addr3) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr4 != -1) {
            if ( calAddrv_addr4( t_Start - prev_t_Start_v_addr4 + prev_t_End_v_addr4, i_Start - prev_i_Start_v_addr4 + prev_i_End_v_addr4, j_Start - prev_j_Start_v_addr4 + prev_j_End_v_addr4) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr4);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr3 = cnt;
                            prev_t_Start_v_addr3 = t_Start;
                            prev_t_End_v_addr3 = t;
                            prev_i_Start_v_addr3 = i_Start;
                            prev_i_End_v_addr3 = i;
                            prev_j_Start_v_addr3 = j_Start;
                            prev_j_End_v_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr4 = cnt;
                            prev_t_Start_v_addr4 = t_Start;
                            prev_t_End_v_addr4 = t;
                            prev_i_Start_v_addr4 = i_Start;
                            prev_i_End_v_addr4 = i;
                            prev_j_Start_v_addr4 = j_Start;
                            prev_j_End_v_addr4 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
            }
            }
            {
            int iLB4 = 1;
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr4(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr5 = -1;
    uint64_t prev_t_Start_v_addr5 = -1;
    uint64_t prev_t_End_v_addr5 = -1;
    uint64_t prev_i_Start_v_addr5 = -1;
    uint64_t prev_i_End_v_addr5 = -1;
    uint64_t prev_j_Start_v_addr5 = -1;
    uint64_t prev_j_End_v_addr5 = -1;
    uint64_t prev_cnt_v_addr6 = -1;
    uint64_t prev_t_Start_v_addr6 = -1;
    uint64_t prev_t_End_v_addr6 = -1;
    uint64_t prev_i_Start_v_addr6 = -1;
    uint64_t prev_i_End_v_addr6 = -1;
    uint64_t prev_j_Start_v_addr6 = -1;
    uint64_t prev_j_End_v_addr6 = -1;
    uint64_t prev_cnt_v_addr7 = -1;
    uint64_t prev_t_Start_v_addr7 = -1;
    uint64_t prev_t_End_v_addr7 = -1;
    uint64_t prev_i_Start_v_addr7 = -1;
    uint64_t prev_i_End_v_addr7 = -1;
    uint64_t prev_j_Start_v_addr7 = -1;
    uint64_t prev_j_End_v_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr5 != -1) {
            if ( calAddrv_addr5( t_Start - prev_t_Start_v_addr5 + prev_t_End_v_addr5, i_Start - prev_i_Start_v_addr5 + prev_i_End_v_addr5, j_Start - prev_j_Start_v_addr5 + prev_j_End_v_addr5) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr6 != -1) {
            if ( calAddrv_addr6( t_Start - prev_t_Start_v_addr6 + prev_t_End_v_addr6, i_Start - prev_i_Start_v_addr6 + prev_i_End_v_addr6, j_Start - prev_j_Start_v_addr6 + prev_j_End_v_addr6) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr7 != -1) {
            if ( calAddrv_addr7( t_Start - prev_t_Start_v_addr7 + prev_t_End_v_addr7, i_Start - prev_i_Start_v_addr7 + prev_i_End_v_addr7, j_Start - prev_j_Start_v_addr7 + prev_j_End_v_addr7) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr5 = cnt;
                            prev_t_Start_v_addr5 = t_Start;
                            prev_t_End_v_addr5 = t;
                            prev_i_Start_v_addr5 = i_Start;
                            prev_i_End_v_addr5 = i;
                            prev_j_Start_v_addr5 = j_Start;
                            prev_j_End_v_addr5 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr6 = cnt;
                            prev_t_Start_v_addr6 = t_Start;
                            prev_t_End_v_addr6 = t;
                            prev_i_Start_v_addr6 = i_Start;
                            prev_i_End_v_addr6 = i;
                            prev_j_Start_v_addr6 = j_Start;
                            prev_j_End_v_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr5(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr7 = cnt;
                            prev_t_Start_v_addr7 = t_Start;
                            prev_t_End_v_addr7 = t;
                            prev_i_Start_v_addr7 = i_Start;
                            prev_i_End_v_addr7 = i;
                            prev_j_Start_v_addr7 = j_Start;
                            prev_j_End_v_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr5 = -1;
    uint64_t prev_t_Start_v_addr5 = -1;
    uint64_t prev_t_End_v_addr5 = -1;
    uint64_t prev_i_Start_v_addr5 = -1;
    uint64_t prev_i_End_v_addr5 = -1;
    uint64_t prev_j_Start_v_addr5 = -1;
    uint64_t prev_j_End_v_addr5 = -1;
    uint64_t prev_cnt_v_addr6 = -1;
    uint64_t prev_t_Start_v_addr6 = -1;
    uint64_t prev_t_End_v_addr6 = -1;
    uint64_t prev_i_Start_v_addr6 = -1;
    uint64_t prev_i_End_v_addr6 = -1;
    uint64_t prev_j_Start_v_addr6 = -1;
    uint64_t prev_j_End_v_addr6 = -1;
    uint64_t prev_cnt_v_addr7 = -1;
    uint64_t prev_t_Start_v_addr7 = -1;
    uint64_t prev_t_End_v_addr7 = -1;
    uint64_t prev_i_Start_v_addr7 = -1;
    uint64_t prev_i_End_v_addr7 = -1;
    uint64_t prev_j_Start_v_addr7 = -1;
    uint64_t prev_j_End_v_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr5 != -1) {
            if ( calAddrv_addr5( t_Start - prev_t_Start_v_addr5 + prev_t_End_v_addr5, i_Start - prev_i_Start_v_addr5 + prev_i_End_v_addr5, j_Start - prev_j_Start_v_addr5 + prev_j_End_v_addr5) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr6 != -1) {
            if ( calAddrv_addr6( t_Start - prev_t_Start_v_addr6 + prev_t_End_v_addr6, i_Start - prev_i_Start_v_addr6 + prev_i_End_v_addr6, j_Start - prev_j_Start_v_addr6 + prev_j_End_v_addr6) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr7 != -1) {
            if ( calAddrv_addr7( t_Start - prev_t_Start_v_addr7 + prev_t_End_v_addr7, i_Start - prev_i_Start_v_addr7 + prev_i_End_v_addr7, j_Start - prev_j_Start_v_addr7 + prev_j_End_v_addr7) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr5 = cnt;
                            prev_t_Start_v_addr5 = t_Start;
                            prev_t_End_v_addr5 = t;
                            prev_i_Start_v_addr5 = i_Start;
                            prev_i_End_v_addr5 = i;
                            prev_j_Start_v_addr5 = j_Start;
                            prev_j_End_v_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr6 = cnt;
                            prev_t_Start_v_addr6 = t_Start;
                            prev_t_End_v_addr6 = t;
                            prev_i_Start_v_addr6 = i_Start;
                            prev_i_End_v_addr6 = i;
                            prev_j_Start_v_addr6 = j_Start;
                            prev_j_End_v_addr6 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr6(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr7 = cnt;
                            prev_t_Start_v_addr7 = t_Start;
                            prev_t_End_v_addr7 = t;
                            prev_i_Start_v_addr7 = i_Start;
                            prev_i_End_v_addr7 = i;
                            prev_j_Start_v_addr7 = j_Start;
                            prev_j_End_v_addr7 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_v_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_v_addr5 = -1;
    uint64_t prev_t_Start_v_addr5 = -1;
    uint64_t prev_t_End_v_addr5 = -1;
    uint64_t prev_i_Start_v_addr5 = -1;
    uint64_t prev_i_End_v_addr5 = -1;
    uint64_t prev_j_Start_v_addr5 = -1;
    uint64_t prev_j_End_v_addr5 = -1;
    uint64_t prev_cnt_v_addr6 = -1;
    uint64_t prev_t_Start_v_addr6 = -1;
    uint64_t prev_t_End_v_addr6 = -1;
    uint64_t prev_i_Start_v_addr6 = -1;
    uint64_t prev_i_End_v_addr6 = -1;
    uint64_t prev_j_Start_v_addr6 = -1;
    uint64_t prev_j_End_v_addr6 = -1;
    uint64_t prev_cnt_v_addr7 = -1;
    uint64_t prev_t_Start_v_addr7 = -1;
    uint64_t prev_t_End_v_addr7 = -1;
    uint64_t prev_i_Start_v_addr7 = -1;
    uint64_t prev_i_End_v_addr7 = -1;
    uint64_t prev_j_Start_v_addr7 = -1;
    uint64_t prev_j_End_v_addr7 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 1) + 1;
        if ( (1023 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_v_addr5 != -1) {
            if ( calAddrv_addr5( t_Start - prev_t_Start_v_addr5 + prev_t_End_v_addr5, i_Start - prev_i_Start_v_addr5 + prev_i_End_v_addr5, j_Start - prev_j_Start_v_addr5 + prev_j_End_v_addr5) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr6 != -1) {
            if ( calAddrv_addr6( t_Start - prev_t_Start_v_addr6 + prev_t_End_v_addr6, i_Start - prev_i_Start_v_addr6 + prev_i_End_v_addr6, j_Start - prev_j_Start_v_addr6 + prev_j_End_v_addr6) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr6);
                goto EndSample;
            }
        }
        if ( prev_cnt_v_addr7 != -1) {
            if ( calAddrv_addr7( t_Start - prev_t_Start_v_addr7 + prev_t_End_v_addr7, i_Start - prev_i_Start_v_addr7 + prev_i_End_v_addr7, j_Start - prev_j_Start_v_addr7 + prev_j_End_v_addr7) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_v_addr7);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int tLB0 = t_Start;
        for ( int t = tLB0; t <= 10; t++) {
            {
            int iLB1 = 1;
            for ( int i = iLB1; i < 1023; i++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr0( t, i) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr1( t, i) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) cnt++;
                {
                int jLB2 = 1;
                for ( int j = jLB2; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv_addr2( t, i) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                {
                int jLB3 = 1022;
                for ( int j = jLB3; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr3( t, i, j) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr4( t, i, j) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
            }
            }
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1023; i++) {
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                if (cntStart == true) cnt++;
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1023; j++) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr5( t, i, j) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr5 = cnt;
                            prev_t_Start_v_addr5 = t_Start;
                            prev_t_End_v_addr5 = t;
                            prev_i_Start_v_addr5 = i_Start;
                            prev_i_End_v_addr5 = i;
                            prev_j_Start_v_addr5 = j_Start;
                            prev_j_End_v_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr6( t, i, j) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr6 = cnt;
                            prev_t_Start_v_addr6 = t_Start;
                            prev_t_End_v_addr6 = t;
                            prev_i_Start_v_addr6 = i_Start;
                            prev_i_End_v_addr6 = i;
                            prev_j_Start_v_addr6 = j_Start;
                            prev_j_End_v_addr6 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrv_addr7( t, i, j) == calAddrv_addr7(t_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_v_addr7 = cnt;
                            prev_t_Start_v_addr7 = t_Start;
                            prev_t_End_v_addr7 = t;
                            prev_i_Start_v_addr7 = i_Start;
                            prev_i_End_v_addr7 = i;
                            prev_j_Start_v_addr7 = j_Start;
                            prev_j_End_v_addr7 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
                if (cntStart == true) cnt++;
                {
                int jLB6 = 1022;
                for ( int j = jLB6; j >= 1; j--) {
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                    if (cntStart == true) cnt++;
                }
                }
            }
            }
        }
        }
EndSample:
        s++;
        }
}
int main() {
    ref_p_addr0();
    ref_p_addr1();
    ref_p_addr2();
    ref_p_addr3();
    ref_p_addr4();
    ref_p_addr5();
    ref_p_addr6();
    ref_p_addr7();
    ref_p_addr8();
    ref_p_addr9();
    ref_q_addr0();
    ref_q_addr1();
    ref_q_addr2();
    ref_q_addr3();
    ref_q_addr4();
    ref_q_addr5();
    ref_q_addr6();
    ref_q_addr7();
    ref_u_addr0();
    ref_u_addr1();
    ref_u_addr2();
    ref_u_addr3();
    ref_u_addr4();
    ref_u_addr5();
    ref_u_addr6();
    ref_u_addr7();
    ref_v_addr0();
    ref_v_addr1();
    ref_v_addr2();
    ref_v_addr3();
    ref_v_addr4();
    ref_v_addr5();
    ref_v_addr6();
    ref_v_addr7();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
adi */ 
