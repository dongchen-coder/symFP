
 /* Start to analysis array index
Array index info
A.addr ((j * 1024) + k)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((j * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + k)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)
A.addr ((i * 1024) + i)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, i)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, j)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access A.addr ((j * 1024) + k)
--------array access A.addr ((i * 1024) + j)
--------array access A.addr ((i * 1024) + j)
------array access A.addr ((j * 1024) + j)
------array access A.addr ((i * 1024) + j)
------array access A.addr ((i * 1024) + j)
----k
----Loop Bound: (0, i)
----Loop inc: (k + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + k)
------array access A.addr ((i * 1024) + i)
------array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)
----array access A.addr ((i * 1024) + i)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
init counter: 0 0 
Dump stride: 1 1 
init counter: 0 0 0 
Dump stride: 1 1 1 
init counter: 0 0 
Dump stride: 1 1 
Dump tree:
----Sample number: 10
------Sample number: 52
--------Sample number: 178
------Sample number: 52
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
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr1( int i, int j, int k) {
    int result = (((j * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr2( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr4( int i, int j) {
    int result = (((j * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr6( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrA_addr7( int i, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr8( int i, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrA_addr9( int i, int k) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrA_addr10( int i, int k) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrA_addr11( int i) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
int calAddrA_addr12( int i) {
    int result = (((i * 1024) + i)) * 8 / 64;
    return result;
}
void ref_A_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_k_Start_A_addr2 = -1;
    uint64_t prev_k_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_k_Start_A_addr3 = -1;
    uint64_t prev_k_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0, k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1, k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2, k_Start - prev_k_Start_A_addr2 + prev_k_End_A_addr2) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3, k_Start - prev_k_Start_A_addr3 + prev_k_End_A_addr3) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            prev_k_Start_A_addr0 = k_Start;
                            prev_k_End_A_addr0 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            prev_k_Start_A_addr1 = k_Start;
                            prev_k_End_A_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            prev_k_Start_A_addr2 = k_Start;
                            prev_k_End_A_addr2 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            prev_k_Start_A_addr3 = k_Start;
                            prev_k_End_A_addr3 = k;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr0(i_Start, j_Start, k_Start)) {
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
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_k_Start_A_addr2 = -1;
    uint64_t prev_k_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_k_Start_A_addr3 = -1;
    uint64_t prev_k_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0, k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1, k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2, k_Start - prev_k_Start_A_addr2 + prev_k_End_A_addr2) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3, k_Start - prev_k_Start_A_addr3 + prev_k_End_A_addr3) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            prev_k_Start_A_addr0 = k_Start;
                            prev_k_End_A_addr0 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            prev_k_Start_A_addr1 = k_Start;
                            prev_k_End_A_addr1 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            prev_k_Start_A_addr2 = k_Start;
                            prev_k_End_A_addr2 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            prev_k_Start_A_addr3 = k_Start;
                            prev_k_End_A_addr3 = k;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr1(i_Start, j_Start, k_Start)) {
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
void ref_A_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_k_Start_A_addr2 = -1;
    uint64_t prev_k_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_k_Start_A_addr3 = -1;
    uint64_t prev_k_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0, k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1, k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2, k_Start - prev_k_Start_A_addr2 + prev_k_End_A_addr2) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3, k_Start - prev_k_Start_A_addr3 + prev_k_End_A_addr3) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            prev_k_Start_A_addr0 = k_Start;
                            prev_k_End_A_addr0 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            prev_k_Start_A_addr1 = k_Start;
                            prev_k_End_A_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            prev_k_Start_A_addr2 = k_Start;
                            prev_k_End_A_addr2 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            prev_k_Start_A_addr3 = k_Start;
                            prev_k_End_A_addr3 = k;
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr2(i_Start, j_Start, k_Start)) {
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
void ref_A_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr0 = -1;
    uint64_t prev_i_Start_A_addr0 = -1;
    uint64_t prev_i_End_A_addr0 = -1;
    uint64_t prev_j_Start_A_addr0 = -1;
    uint64_t prev_j_End_A_addr0 = -1;
    uint64_t prev_k_Start_A_addr0 = -1;
    uint64_t prev_k_End_A_addr0 = -1;
    uint64_t prev_cnt_A_addr1 = -1;
    uint64_t prev_i_Start_A_addr1 = -1;
    uint64_t prev_i_End_A_addr1 = -1;
    uint64_t prev_j_Start_A_addr1 = -1;
    uint64_t prev_j_End_A_addr1 = -1;
    uint64_t prev_k_Start_A_addr1 = -1;
    uint64_t prev_k_End_A_addr1 = -1;
    uint64_t prev_cnt_A_addr2 = -1;
    uint64_t prev_i_Start_A_addr2 = -1;
    uint64_t prev_i_End_A_addr2 = -1;
    uint64_t prev_j_Start_A_addr2 = -1;
    uint64_t prev_j_End_A_addr2 = -1;
    uint64_t prev_k_Start_A_addr2 = -1;
    uint64_t prev_k_End_A_addr2 = -1;
    uint64_t prev_cnt_A_addr3 = -1;
    uint64_t prev_i_Start_A_addr3 = -1;
    uint64_t prev_i_End_A_addr3 = -1;
    uint64_t prev_j_Start_A_addr3 = -1;
    uint64_t prev_j_End_A_addr3 = -1;
    uint64_t prev_k_Start_A_addr3 = -1;
    uint64_t prev_k_End_A_addr3 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 178;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        if ( (j_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (j_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr0 != -1) {
            if ( calAddrA_addr0( i_Start - prev_i_Start_A_addr0 + prev_i_End_A_addr0, j_Start - prev_j_Start_A_addr0 + prev_j_End_A_addr0, k_Start - prev_k_Start_A_addr0 + prev_k_End_A_addr0) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr1 != -1) {
            if ( calAddrA_addr1( i_Start - prev_i_Start_A_addr1 + prev_i_End_A_addr1, j_Start - prev_j_Start_A_addr1 + prev_j_End_A_addr1, k_Start - prev_k_Start_A_addr1 + prev_k_End_A_addr1) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr2 != -1) {
            if ( calAddrA_addr2( i_Start - prev_i_Start_A_addr2 + prev_i_End_A_addr2, j_Start - prev_j_Start_A_addr2 + prev_j_End_A_addr2, k_Start - prev_k_Start_A_addr2 + prev_k_End_A_addr2) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr3 != -1) {
            if ( calAddrA_addr3( i_Start - prev_i_Start_A_addr3 + prev_i_End_A_addr3, j_Start - prev_j_Start_A_addr3 + prev_j_End_A_addr3, k_Start - prev_k_Start_A_addr3 + prev_k_End_A_addr3) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr3);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                if ( i == i_Start && j == j_Start ) {
                    kLB2 = k_Start;
                }
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr0 = cnt;
                            prev_i_Start_A_addr0 = i_Start;
                            prev_i_End_A_addr0 = i;
                            prev_j_Start_A_addr0 = j_Start;
                            prev_j_End_A_addr0 = j;
                            prev_k_Start_A_addr0 = k_Start;
                            prev_k_End_A_addr0 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr1 = cnt;
                            prev_i_Start_A_addr1 = i_Start;
                            prev_i_End_A_addr1 = i;
                            prev_j_Start_A_addr1 = j_Start;
                            prev_j_End_A_addr1 = j;
                            prev_k_Start_A_addr1 = k_Start;
                            prev_k_End_A_addr1 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr2 = cnt;
                            prev_i_Start_A_addr2 = i_Start;
                            prev_i_End_A_addr2 = i;
                            prev_j_Start_A_addr2 = j_Start;
                            prev_j_End_A_addr2 = j;
                            prev_k_Start_A_addr2 = k_Start;
                            prev_k_End_A_addr2 = k;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_A_addr3 = cnt;
                            prev_i_Start_A_addr3 = i_Start;
                            prev_i_End_A_addr3 = i;
                            prev_j_Start_A_addr3 = j_Start;
                            prev_j_End_A_addr3 = j;
                            prev_k_Start_A_addr3 = k_Start;
                            prev_k_End_A_addr3 = k;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr3(i_Start, j_Start, k_Start)) {
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
void ref_A_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    uint64_t prev_cnt_A_addr6 = -1;
    uint64_t prev_i_Start_A_addr6 = -1;
    uint64_t prev_i_End_A_addr6 = -1;
    uint64_t prev_j_Start_A_addr6 = -1;
    uint64_t prev_j_End_A_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5) == calAddrA_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr6 != -1) {
            if ( calAddrA_addr6( i_Start - prev_i_Start_A_addr6 + prev_i_End_A_addr6, j_Start - prev_j_Start_A_addr6 + prev_j_End_A_addr6) == calAddrA_addr4(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr4(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr4 = cnt;
                        prev_i_Start_A_addr4 = i_Start;
                        prev_i_End_A_addr4 = i;
                        prev_j_Start_A_addr4 = j_Start;
                        prev_j_End_A_addr4 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr5 = cnt;
                        prev_i_Start_A_addr5 = i_Start;
                        prev_i_End_A_addr5 = i;
                        prev_j_Start_A_addr5 = j_Start;
                        prev_j_End_A_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr6 = cnt;
                        prev_i_Start_A_addr6 = i_Start;
                        prev_i_End_A_addr6 = i;
                        prev_j_Start_A_addr6 = j_Start;
                        prev_j_End_A_addr6 = j;
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr4(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr4(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr4(i_Start, j_Start)) {
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
void ref_A_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    uint64_t prev_cnt_A_addr6 = -1;
    uint64_t prev_i_Start_A_addr6 = -1;
    uint64_t prev_i_End_A_addr6 = -1;
    uint64_t prev_j_Start_A_addr6 = -1;
    uint64_t prev_j_End_A_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5) == calAddrA_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr6 != -1) {
            if ( calAddrA_addr6( i_Start - prev_i_Start_A_addr6 + prev_i_End_A_addr6, j_Start - prev_j_Start_A_addr6 + prev_j_End_A_addr6) == calAddrA_addr5(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr5(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr4 = cnt;
                        prev_i_Start_A_addr4 = i_Start;
                        prev_i_End_A_addr4 = i;
                        prev_j_Start_A_addr4 = j_Start;
                        prev_j_End_A_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr5 = cnt;
                        prev_i_Start_A_addr5 = i_Start;
                        prev_i_End_A_addr5 = i;
                        prev_j_Start_A_addr5 = j_Start;
                        prev_j_End_A_addr5 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr6 = cnt;
                        prev_i_Start_A_addr6 = i_Start;
                        prev_i_End_A_addr6 = i;
                        prev_j_Start_A_addr6 = j_Start;
                        prev_j_End_A_addr6 = j;
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr5(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr5(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr5(i_Start, j_Start)) {
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
void ref_A_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr4 = -1;
    uint64_t prev_i_Start_A_addr4 = -1;
    uint64_t prev_i_End_A_addr4 = -1;
    uint64_t prev_j_Start_A_addr4 = -1;
    uint64_t prev_j_End_A_addr4 = -1;
    uint64_t prev_cnt_A_addr5 = -1;
    uint64_t prev_i_Start_A_addr5 = -1;
    uint64_t prev_i_End_A_addr5 = -1;
    uint64_t prev_j_Start_A_addr5 = -1;
    uint64_t prev_j_End_A_addr5 = -1;
    uint64_t prev_cnt_A_addr6 = -1;
    uint64_t prev_i_Start_A_addr6 = -1;
    uint64_t prev_i_End_A_addr6 = -1;
    uint64_t prev_j_Start_A_addr6 = -1;
    uint64_t prev_j_End_A_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr4 != -1) {
            if ( calAddrA_addr4( i_Start - prev_i_Start_A_addr4 + prev_i_End_A_addr4, j_Start - prev_j_Start_A_addr4 + prev_j_End_A_addr4) == calAddrA_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr5 != -1) {
            if ( calAddrA_addr5( i_Start - prev_i_Start_A_addr5 + prev_i_End_A_addr5, j_Start - prev_j_Start_A_addr5 + prev_j_End_A_addr5) == calAddrA_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr6 != -1) {
            if ( calAddrA_addr6( i_Start - prev_i_Start_A_addr6 + prev_i_End_A_addr6, j_Start - prev_j_Start_A_addr6 + prev_j_End_A_addr6) == calAddrA_addr6(i_Start, j_Start)) {
                rtHistoCal(prev_cnt_A_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr6(i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr4 = cnt;
                        prev_i_Start_A_addr4 = i_Start;
                        prev_i_End_A_addr4 = i;
                        prev_j_Start_A_addr4 = j_Start;
                        prev_j_End_A_addr4 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr5 = cnt;
                        prev_i_Start_A_addr5 = i_Start;
                        prev_i_End_A_addr5 = i;
                        prev_j_Start_A_addr5 = j_Start;
                        prev_j_End_A_addr5 = j;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr6 = cnt;
                        prev_i_Start_A_addr6 = i_Start;
                        prev_i_End_A_addr6 = i;
                        prev_j_Start_A_addr6 = j_Start;
                        prev_j_End_A_addr6 = j;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr6(i_Start, j_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr6(i_Start, j_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr6(i_Start, j_Start)) {
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
void ref_A_addr7() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr7 = -1;
    uint64_t prev_i_Start_A_addr7 = -1;
    uint64_t prev_i_End_A_addr7 = -1;
    uint64_t prev_k_Start_A_addr7 = -1;
    uint64_t prev_k_End_A_addr7 = -1;
    uint64_t prev_cnt_A_addr8 = -1;
    uint64_t prev_i_Start_A_addr8 = -1;
    uint64_t prev_i_End_A_addr8 = -1;
    uint64_t prev_k_Start_A_addr8 = -1;
    uint64_t prev_k_End_A_addr8 = -1;
    uint64_t prev_cnt_A_addr9 = -1;
    uint64_t prev_i_Start_A_addr9 = -1;
    uint64_t prev_i_End_A_addr9 = -1;
    uint64_t prev_k_Start_A_addr9 = -1;
    uint64_t prev_k_End_A_addr9 = -1;
    uint64_t prev_cnt_A_addr10 = -1;
    uint64_t prev_i_Start_A_addr10 = -1;
    uint64_t prev_i_End_A_addr10 = -1;
    uint64_t prev_k_Start_A_addr10 = -1;
    uint64_t prev_k_End_A_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr7 != -1) {
            if ( calAddrA_addr7( i_Start - prev_i_Start_A_addr7 + prev_i_End_A_addr7, k_Start - prev_k_Start_A_addr7 + prev_k_End_A_addr7) == calAddrA_addr7(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr8 != -1) {
            if ( calAddrA_addr8( i_Start - prev_i_Start_A_addr8 + prev_i_End_A_addr8, k_Start - prev_k_Start_A_addr8 + prev_k_End_A_addr8) == calAddrA_addr7(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr9 != -1) {
            if ( calAddrA_addr9( i_Start - prev_i_Start_A_addr9 + prev_i_End_A_addr9, k_Start - prev_k_Start_A_addr9 + prev_k_End_A_addr9) == calAddrA_addr7(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr10 != -1) {
            if ( calAddrA_addr10( i_Start - prev_i_Start_A_addr10 + prev_i_End_A_addr10, k_Start - prev_k_Start_A_addr10 + prev_k_End_A_addr10) == calAddrA_addr7(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr7(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr7 = cnt;
                        prev_i_Start_A_addr7 = i_Start;
                        prev_i_End_A_addr7 = i;
                        prev_k_Start_A_addr7 = k_Start;
                        prev_k_End_A_addr7 = k;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr8 = cnt;
                        prev_i_Start_A_addr8 = i_Start;
                        prev_i_End_A_addr8 = i;
                        prev_k_Start_A_addr8 = k_Start;
                        prev_k_End_A_addr8 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr9 = cnt;
                        prev_i_Start_A_addr9 = i_Start;
                        prev_i_End_A_addr9 = i;
                        prev_k_Start_A_addr9 = k_Start;
                        prev_k_End_A_addr9 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr7(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr10 = cnt;
                        prev_i_Start_A_addr10 = i_Start;
                        prev_i_End_A_addr10 = i;
                        prev_k_Start_A_addr10 = k_Start;
                        prev_k_End_A_addr10 = k;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr7(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr7(i_Start, k_Start)) {
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
void ref_A_addr8() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr7 = -1;
    uint64_t prev_i_Start_A_addr7 = -1;
    uint64_t prev_i_End_A_addr7 = -1;
    uint64_t prev_k_Start_A_addr7 = -1;
    uint64_t prev_k_End_A_addr7 = -1;
    uint64_t prev_cnt_A_addr8 = -1;
    uint64_t prev_i_Start_A_addr8 = -1;
    uint64_t prev_i_End_A_addr8 = -1;
    uint64_t prev_k_Start_A_addr8 = -1;
    uint64_t prev_k_End_A_addr8 = -1;
    uint64_t prev_cnt_A_addr9 = -1;
    uint64_t prev_i_Start_A_addr9 = -1;
    uint64_t prev_i_End_A_addr9 = -1;
    uint64_t prev_k_Start_A_addr9 = -1;
    uint64_t prev_k_End_A_addr9 = -1;
    uint64_t prev_cnt_A_addr10 = -1;
    uint64_t prev_i_Start_A_addr10 = -1;
    uint64_t prev_i_End_A_addr10 = -1;
    uint64_t prev_k_Start_A_addr10 = -1;
    uint64_t prev_k_End_A_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr7 != -1) {
            if ( calAddrA_addr7( i_Start - prev_i_Start_A_addr7 + prev_i_End_A_addr7, k_Start - prev_k_Start_A_addr7 + prev_k_End_A_addr7) == calAddrA_addr8(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr8 != -1) {
            if ( calAddrA_addr8( i_Start - prev_i_Start_A_addr8 + prev_i_End_A_addr8, k_Start - prev_k_Start_A_addr8 + prev_k_End_A_addr8) == calAddrA_addr8(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr9 != -1) {
            if ( calAddrA_addr9( i_Start - prev_i_Start_A_addr9 + prev_i_End_A_addr9, k_Start - prev_k_Start_A_addr9 + prev_k_End_A_addr9) == calAddrA_addr8(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr10 != -1) {
            if ( calAddrA_addr10( i_Start - prev_i_Start_A_addr10 + prev_i_End_A_addr10, k_Start - prev_k_Start_A_addr10 + prev_k_End_A_addr10) == calAddrA_addr8(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr8(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr7 = cnt;
                        prev_i_Start_A_addr7 = i_Start;
                        prev_i_End_A_addr7 = i;
                        prev_k_Start_A_addr7 = k_Start;
                        prev_k_End_A_addr7 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr8 = cnt;
                        prev_i_Start_A_addr8 = i_Start;
                        prev_i_End_A_addr8 = i;
                        prev_k_Start_A_addr8 = k_Start;
                        prev_k_End_A_addr8 = k;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr9 = cnt;
                        prev_i_Start_A_addr9 = i_Start;
                        prev_i_End_A_addr9 = i;
                        prev_k_Start_A_addr9 = k_Start;
                        prev_k_End_A_addr9 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr8(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr10 = cnt;
                        prev_i_Start_A_addr10 = i_Start;
                        prev_i_End_A_addr10 = i;
                        prev_k_Start_A_addr10 = k_Start;
                        prev_k_End_A_addr10 = k;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr8(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr8(i_Start, k_Start)) {
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
void ref_A_addr9() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr7 = -1;
    uint64_t prev_i_Start_A_addr7 = -1;
    uint64_t prev_i_End_A_addr7 = -1;
    uint64_t prev_k_Start_A_addr7 = -1;
    uint64_t prev_k_End_A_addr7 = -1;
    uint64_t prev_cnt_A_addr8 = -1;
    uint64_t prev_i_Start_A_addr8 = -1;
    uint64_t prev_i_End_A_addr8 = -1;
    uint64_t prev_k_Start_A_addr8 = -1;
    uint64_t prev_k_End_A_addr8 = -1;
    uint64_t prev_cnt_A_addr9 = -1;
    uint64_t prev_i_Start_A_addr9 = -1;
    uint64_t prev_i_End_A_addr9 = -1;
    uint64_t prev_k_Start_A_addr9 = -1;
    uint64_t prev_k_End_A_addr9 = -1;
    uint64_t prev_cnt_A_addr10 = -1;
    uint64_t prev_i_Start_A_addr10 = -1;
    uint64_t prev_i_End_A_addr10 = -1;
    uint64_t prev_k_Start_A_addr10 = -1;
    uint64_t prev_k_End_A_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr7 != -1) {
            if ( calAddrA_addr7( i_Start - prev_i_Start_A_addr7 + prev_i_End_A_addr7, k_Start - prev_k_Start_A_addr7 + prev_k_End_A_addr7) == calAddrA_addr9(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr8 != -1) {
            if ( calAddrA_addr8( i_Start - prev_i_Start_A_addr8 + prev_i_End_A_addr8, k_Start - prev_k_Start_A_addr8 + prev_k_End_A_addr8) == calAddrA_addr9(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr9 != -1) {
            if ( calAddrA_addr9( i_Start - prev_i_Start_A_addr9 + prev_i_End_A_addr9, k_Start - prev_k_Start_A_addr9 + prev_k_End_A_addr9) == calAddrA_addr9(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr10 != -1) {
            if ( calAddrA_addr10( i_Start - prev_i_Start_A_addr10 + prev_i_End_A_addr10, k_Start - prev_k_Start_A_addr10 + prev_k_End_A_addr10) == calAddrA_addr9(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr9(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr7 = cnt;
                        prev_i_Start_A_addr7 = i_Start;
                        prev_i_End_A_addr7 = i;
                        prev_k_Start_A_addr7 = k_Start;
                        prev_k_End_A_addr7 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr8 = cnt;
                        prev_i_Start_A_addr8 = i_Start;
                        prev_i_End_A_addr8 = i;
                        prev_k_Start_A_addr8 = k_Start;
                        prev_k_End_A_addr8 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr9 = cnt;
                        prev_i_Start_A_addr9 = i_Start;
                        prev_i_End_A_addr9 = i;
                        prev_k_Start_A_addr9 = k_Start;
                        prev_k_End_A_addr9 = k;
                        goto EndSample;
                    }
                }
                cntStart = true;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr9(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr10 = cnt;
                        prev_i_Start_A_addr10 = i_Start;
                        prev_i_End_A_addr10 = i;
                        prev_k_Start_A_addr10 = k_Start;
                        prev_k_End_A_addr10 = k;
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr9(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr9(i_Start, k_Start)) {
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
void ref_A_addr10() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr7 = -1;
    uint64_t prev_i_Start_A_addr7 = -1;
    uint64_t prev_i_End_A_addr7 = -1;
    uint64_t prev_k_Start_A_addr7 = -1;
    uint64_t prev_k_End_A_addr7 = -1;
    uint64_t prev_cnt_A_addr8 = -1;
    uint64_t prev_i_Start_A_addr8 = -1;
    uint64_t prev_i_End_A_addr8 = -1;
    uint64_t prev_k_Start_A_addr8 = -1;
    uint64_t prev_k_End_A_addr8 = -1;
    uint64_t prev_cnt_A_addr9 = -1;
    uint64_t prev_i_Start_A_addr9 = -1;
    uint64_t prev_i_End_A_addr9 = -1;
    uint64_t prev_k_Start_A_addr9 = -1;
    uint64_t prev_k_End_A_addr9 = -1;
    uint64_t prev_cnt_A_addr10 = -1;
    uint64_t prev_i_Start_A_addr10 = -1;
    uint64_t prev_i_End_A_addr10 = -1;
    uint64_t prev_k_Start_A_addr10 = -1;
    uint64_t prev_k_End_A_addr10 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 52;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (i_Start - 0) == 0) goto SAMPLE;
        int k_Start = rand() % (i_Start - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr7 != -1) {
            if ( calAddrA_addr7( i_Start - prev_i_Start_A_addr7 + prev_i_End_A_addr7, k_Start - prev_k_Start_A_addr7 + prev_k_End_A_addr7) == calAddrA_addr10(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr7);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr8 != -1) {
            if ( calAddrA_addr8( i_Start - prev_i_Start_A_addr8 + prev_i_End_A_addr8, k_Start - prev_k_Start_A_addr8 + prev_k_End_A_addr8) == calAddrA_addr10(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr8);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr9 != -1) {
            if ( calAddrA_addr9( i_Start - prev_i_Start_A_addr9 + prev_i_End_A_addr9, k_Start - prev_k_Start_A_addr9 + prev_k_End_A_addr9) == calAddrA_addr10(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr9);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr10 != -1) {
            if ( calAddrA_addr10( i_Start - prev_i_Start_A_addr10 + prev_i_End_A_addr10, k_Start - prev_k_Start_A_addr10 + prev_k_End_A_addr10) == calAddrA_addr10(i_Start, k_Start)) {
                rtHistoCal(prev_cnt_A_addr10);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr10(i_Start, k_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            if ( i == i_Start ) {
                kLB3 = k_Start;
            }
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr7 = cnt;
                        prev_i_Start_A_addr7 = i_Start;
                        prev_i_End_A_addr7 = i;
                        prev_k_Start_A_addr7 = k_Start;
                        prev_k_End_A_addr7 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr8 = cnt;
                        prev_i_Start_A_addr8 = i_Start;
                        prev_i_End_A_addr8 = i;
                        prev_k_Start_A_addr8 = k_Start;
                        prev_k_End_A_addr8 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr9 = cnt;
                        prev_i_Start_A_addr9 = i_Start;
                        prev_i_End_A_addr9 = i;
                        prev_k_Start_A_addr9 = k_Start;
                        prev_k_End_A_addr9 = k;
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr10(i_Start, k_Start)) {
                        rtHistoCal(cnt);
                        prev_cnt_A_addr10 = cnt;
                        prev_i_Start_A_addr10 = i_Start;
                        prev_i_End_A_addr10 = i;
                        prev_k_Start_A_addr10 = k_Start;
                        prev_k_End_A_addr10 = k;
                        goto EndSample;
                    }
                }
                cntStart = true;
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr10(i_Start, k_Start)) {
                    rtHistoCal(cnt);
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr10(i_Start, k_Start)) {
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
void ref_A_addr11() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr11 = -1;
    uint64_t prev_i_Start_A_addr11 = -1;
    uint64_t prev_i_End_A_addr11 = -1;
    uint64_t prev_cnt_A_addr12 = -1;
    uint64_t prev_i_Start_A_addr12 = -1;
    uint64_t prev_i_End_A_addr12 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr11 != -1) {
            if ( calAddrA_addr11( i_Start - prev_i_Start_A_addr11 + prev_i_End_A_addr11) == calAddrA_addr11(i_Start)) {
                rtHistoCal(prev_cnt_A_addr11);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr12 != -1) {
            if ( calAddrA_addr12( i_Start - prev_i_Start_A_addr12 + prev_i_End_A_addr12) == calAddrA_addr11(i_Start)) {
                rtHistoCal(prev_cnt_A_addr12);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr11(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr11(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr11(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_A_addr11 = cnt;
                    prev_i_Start_A_addr11 = i_Start;
                    prev_i_End_A_addr11 = i;
                    goto EndSample;
                }
            }
            cntStart = true;
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr11(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_A_addr12 = cnt;
                    prev_i_Start_A_addr12 = i_Start;
                    prev_i_End_A_addr12 = i;
                    goto EndSample;
                }
            }
        }
        }
EndSample:
        s++;
        }
}
void ref_A_addr12() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_A_addr11 = -1;
    uint64_t prev_i_Start_A_addr11 = -1;
    uint64_t prev_i_End_A_addr11 = -1;
    uint64_t prev_cnt_A_addr12 = -1;
    uint64_t prev_i_Start_A_addr12 = -1;
    uint64_t prev_i_End_A_addr12 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_A_addr11 != -1) {
            if ( calAddrA_addr11( i_Start - prev_i_Start_A_addr11 + prev_i_End_A_addr11) == calAddrA_addr12(i_Start)) {
                rtHistoCal(prev_cnt_A_addr11);
                goto EndSample;
            }
        }
        if ( prev_cnt_A_addr12 != -1) {
            if ( calAddrA_addr12( i_Start - prev_i_Start_A_addr12 + prev_i_End_A_addr12) == calAddrA_addr12(i_Start)) {
                rtHistoCal(prev_cnt_A_addr12);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 1024; i++) {
            {
            int jLB1 = 0;
            for ( int j = jLB1; j < i; j++) {
                {
                int kLB2 = 0;
                for ( int k = kLB2; k < j; k++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr2( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr3( i, j, k) == calAddrA_addr12(i_Start)) {
                            rtHistoCal(cnt);
                            goto EndSample;
                        }
                    }
                }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr4( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr5( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr6( i, j) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            {
            int kLB3 = 0;
            for ( int k = kLB3; k < i; k++) {
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr7( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr8( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr9( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr10( i, k) == calAddrA_addr12(i_Start)) {
                        rtHistoCal(cnt);
                        goto EndSample;
                    }
                }
            }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr11( i) == calAddrA_addr12(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_A_addr11 = cnt;
                    prev_i_Start_A_addr11 = i_Start;
                    prev_i_End_A_addr11 = i;
                    goto EndSample;
                }
            }
            if (cntStart == true) {
                cnt++;
                if ( calAddrA_addr12( i) == calAddrA_addr12(i_Start)) {
                    rtHistoCal(cnt);
                    prev_cnt_A_addr12 = cnt;
                    prev_i_Start_A_addr12 = i_Start;
                    prev_i_End_A_addr12 = i;
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
int main() {
    ref_A_addr0();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr3();
    ref_A_addr4();
    ref_A_addr5();
    ref_A_addr6();
    ref_A_addr7();
    ref_A_addr8();
    ref_A_addr9();
    ref_A_addr10();
    ref_A_addr11();
    ref_A_addr12();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
cholesky */ 
