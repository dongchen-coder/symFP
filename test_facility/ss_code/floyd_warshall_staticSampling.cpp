
 /* Start to analysis array index
Array index info
path.addr ((k * 1024) + j)
path.addr ((i * 1024) + j)
path.addr ((i * 1024) + k)
path.addr ((i * 1024) + j)
path.addr ((i * 1024) + k)
path.addr ((k * 1024) + j)
path.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %path

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--k
--Loop Bound: (0, 1024)
--Loop inc: (k + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 1024)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access path.addr ((i * 1024) + j)
--------array access path.addr ((i * 1024) + k)
--------array access path.addr ((k * 1024) + j)
--------array access path.addr ((i * 1024) + j)
--------array access path.addr ((i * 1024) + k)
--------array access path.addr ((k * 1024) + j)
--------array access path.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 10
------Sample number: 104
--------Sample number: 1073
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
int calAddrpath_addr0( int k, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrpath_addr1( int k, int i, int j) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrpath_addr2( int k, int i, int j) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrpath_addr3( int k, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrpath_addr4( int k, int i, int j) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
int calAddrpath_addr5( int k, int i, int j) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
int calAddrpath_addr6( int k, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_path_addr0() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr0(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr1() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr1(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr2() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr2(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr3() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr3(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr4() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr4(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr5() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    cntStart = true;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr5(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
void ref_path_addr6() {
    /* Generating search reuse init code (same loop) */
    uint64_t prev_cnt_path_addr0 = -1;
    uint64_t prev_k_Start_path_addr0 = -1;
    uint64_t prev_k_End_path_addr0 = -1;
    uint64_t prev_i_Start_path_addr0 = -1;
    uint64_t prev_i_End_path_addr0 = -1;
    uint64_t prev_j_Start_path_addr0 = -1;
    uint64_t prev_j_End_path_addr0 = -1;
    uint64_t prev_cnt_path_addr1 = -1;
    uint64_t prev_k_Start_path_addr1 = -1;
    uint64_t prev_k_End_path_addr1 = -1;
    uint64_t prev_i_Start_path_addr1 = -1;
    uint64_t prev_i_End_path_addr1 = -1;
    uint64_t prev_j_Start_path_addr1 = -1;
    uint64_t prev_j_End_path_addr1 = -1;
    uint64_t prev_cnt_path_addr2 = -1;
    uint64_t prev_k_Start_path_addr2 = -1;
    uint64_t prev_k_End_path_addr2 = -1;
    uint64_t prev_i_Start_path_addr2 = -1;
    uint64_t prev_i_End_path_addr2 = -1;
    uint64_t prev_j_Start_path_addr2 = -1;
    uint64_t prev_j_End_path_addr2 = -1;
    uint64_t prev_cnt_path_addr3 = -1;
    uint64_t prev_k_Start_path_addr3 = -1;
    uint64_t prev_k_End_path_addr3 = -1;
    uint64_t prev_i_Start_path_addr3 = -1;
    uint64_t prev_i_End_path_addr3 = -1;
    uint64_t prev_j_Start_path_addr3 = -1;
    uint64_t prev_j_End_path_addr3 = -1;
    uint64_t prev_cnt_path_addr4 = -1;
    uint64_t prev_k_Start_path_addr4 = -1;
    uint64_t prev_k_End_path_addr4 = -1;
    uint64_t prev_i_Start_path_addr4 = -1;
    uint64_t prev_i_End_path_addr4 = -1;
    uint64_t prev_j_Start_path_addr4 = -1;
    uint64_t prev_j_End_path_addr4 = -1;
    uint64_t prev_cnt_path_addr5 = -1;
    uint64_t prev_k_Start_path_addr5 = -1;
    uint64_t prev_k_End_path_addr5 = -1;
    uint64_t prev_i_Start_path_addr5 = -1;
    uint64_t prev_i_End_path_addr5 = -1;
    uint64_t prev_j_Start_path_addr5 = -1;
    uint64_t prev_j_End_path_addr5 = -1;
    uint64_t prev_cnt_path_addr6 = -1;
    uint64_t prev_k_Start_path_addr6 = -1;
    uint64_t prev_k_End_path_addr6 = -1;
    uint64_t prev_i_Start_path_addr6 = -1;
    uint64_t prev_i_End_path_addr6 = -1;
    uint64_t prev_j_Start_path_addr6 = -1;
    uint64_t prev_j_End_path_addr6 = -1;
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1073;) {
SAMPLE:
        int k_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(k_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating search reuse body code (use reuse are in the same loop) */
        if ( prev_cnt_path_addr0 != -1) {
            if ( calAddrpath_addr0( k_Start - prev_k_Start_path_addr0 + prev_k_End_path_addr0, i_Start - prev_i_Start_path_addr0 + prev_i_End_path_addr0, j_Start - prev_j_Start_path_addr0 + prev_j_End_path_addr0) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr0);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr1 != -1) {
            if ( calAddrpath_addr1( k_Start - prev_k_Start_path_addr1 + prev_k_End_path_addr1, i_Start - prev_i_Start_path_addr1 + prev_i_End_path_addr1, j_Start - prev_j_Start_path_addr1 + prev_j_End_path_addr1) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr1);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr2 != -1) {
            if ( calAddrpath_addr2( k_Start - prev_k_Start_path_addr2 + prev_k_End_path_addr2, i_Start - prev_i_Start_path_addr2 + prev_i_End_path_addr2, j_Start - prev_j_Start_path_addr2 + prev_j_End_path_addr2) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr2);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr3 != -1) {
            if ( calAddrpath_addr3( k_Start - prev_k_Start_path_addr3 + prev_k_End_path_addr3, i_Start - prev_i_Start_path_addr3 + prev_i_End_path_addr3, j_Start - prev_j_Start_path_addr3 + prev_j_End_path_addr3) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr3);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr4 != -1) {
            if ( calAddrpath_addr4( k_Start - prev_k_Start_path_addr4 + prev_k_End_path_addr4, i_Start - prev_i_Start_path_addr4 + prev_i_End_path_addr4, j_Start - prev_j_Start_path_addr4 + prev_j_End_path_addr4) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr4);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr5 != -1) {
            if ( calAddrpath_addr5( k_Start - prev_k_Start_path_addr5 + prev_k_End_path_addr5, i_Start - prev_i_Start_path_addr5 + prev_i_End_path_addr5, j_Start - prev_j_Start_path_addr5 + prev_j_End_path_addr5) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr5);
                goto EndSample;
            }
        }
        if ( prev_cnt_path_addr6 != -1) {
            if ( calAddrpath_addr6( k_Start - prev_k_Start_path_addr6 + prev_k_End_path_addr6, i_Start - prev_i_Start_path_addr6 + prev_i_End_path_addr6, j_Start - prev_j_Start_path_addr6 + prev_j_End_path_addr6) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                rtHistoCal(prev_cnt_path_addr6);
                goto EndSample;
            }
        }
        /* Generating reuse search code */

        {
        int kLB0 = k_Start;
        for ( int k = kLB0; k < 1024; k++) {
            {
            int iLB1 = 0;
            if ( k == k_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 1024; i++) {
                {
                int jLB2 = 0;
                if ( k == k_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 1024; j++) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr0( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr0 = cnt;
                            prev_k_Start_path_addr0 = k_Start;
                            prev_k_End_path_addr0 = k;
                            prev_i_Start_path_addr0 = i_Start;
                            prev_i_End_path_addr0 = i;
                            prev_j_Start_path_addr0 = j_Start;
                            prev_j_End_path_addr0 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr1( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr1 = cnt;
                            prev_k_Start_path_addr1 = k_Start;
                            prev_k_End_path_addr1 = k;
                            prev_i_Start_path_addr1 = i_Start;
                            prev_i_End_path_addr1 = i;
                            prev_j_Start_path_addr1 = j_Start;
                            prev_j_End_path_addr1 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr2( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr2 = cnt;
                            prev_k_Start_path_addr2 = k_Start;
                            prev_k_End_path_addr2 = k;
                            prev_i_Start_path_addr2 = i_Start;
                            prev_i_End_path_addr2 = i;
                            prev_j_Start_path_addr2 = j_Start;
                            prev_j_End_path_addr2 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr3( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr3 = cnt;
                            prev_k_Start_path_addr3 = k_Start;
                            prev_k_End_path_addr3 = k;
                            prev_i_Start_path_addr3 = i_Start;
                            prev_i_End_path_addr3 = i;
                            prev_j_Start_path_addr3 = j_Start;
                            prev_j_End_path_addr3 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr4( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr4 = cnt;
                            prev_k_Start_path_addr4 = k_Start;
                            prev_k_End_path_addr4 = k;
                            prev_i_Start_path_addr4 = i_Start;
                            prev_i_End_path_addr4 = i;
                            prev_j_Start_path_addr4 = j_Start;
                            prev_j_End_path_addr4 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr5( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr5 = cnt;
                            prev_k_Start_path_addr5 = k_Start;
                            prev_k_End_path_addr5 = k;
                            prev_i_Start_path_addr5 = i_Start;
                            prev_i_End_path_addr5 = i;
                            prev_j_Start_path_addr5 = j_Start;
                            prev_j_End_path_addr5 = j;
                            goto EndSample;
                        }
                    }
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrpath_addr6( k, i, j) == calAddrpath_addr6(k_Start, i_Start, j_Start)) {
                            rtHistoCal(cnt);
                            prev_cnt_path_addr6 = cnt;
                            prev_k_Start_path_addr6 = k_Start;
                            prev_k_End_path_addr6 = k;
                            prev_i_Start_path_addr6 = i_Start;
                            prev_i_End_path_addr6 = i;
                            prev_j_Start_path_addr6 = j_Start;
                            prev_j_End_path_addr6 = j;
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
int main() {
    ref_path_addr0();
    ref_path_addr1();
    ref_path_addr2();
    ref_path_addr3();
    ref_path_addr4();
    ref_path_addr5();
    ref_path_addr6();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Start to analyze function:  
floyd_warshall */ 
