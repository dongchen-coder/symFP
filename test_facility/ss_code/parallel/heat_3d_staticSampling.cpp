
 /* Start to analysis array index
Array index info: Total number of references: 22
A.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
A.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
A.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
B.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
B.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
A.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
B.addr ((((i * 256) * 256) + (j * 256)) + k)
A.addr ((((i * 256) * 256) + (j * 256)) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %B
double* %A

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (1, 10)
--Loop inc: (t + 1)
--Loop predicate: <=
----i
----Loop Bound: (1, 255)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 255)
------Loop inc: (j + 1)
------Loop predicate: <
--------k
--------Loop Bound: (1, 255)
--------Loop inc: (k + 1)
--------Loop predicate: <
----------array access A.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
----------array access A.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----i
----Loop Bound: (1, 255)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 255)
------Loop inc: (j + 1)
------Loop predicate: <
--------k
--------Loop Bound: (1, 255)
--------Loop inc: (k + 1)
--------Loop predicate: <
----------array access B.addr (((((i + 1) * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr (((((i - 1) * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + ((j + 1) * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr ((((i * 256) * 256) + ((j - 1) * 256)) + k)
----------array access B.addr (((((i * 256) * 256) + (j * 256)) + k) + 1)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access B.addr (((((i * 256) * 256) + (j * 256)) + k) - 1)
----------array access B.addr ((((i * 256) * 256) + (j * 256)) + k)
----------array access A.addr ((((i * 256) * 256) + (j * 256)) + k)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 0
--------Sample number: 0
----------Sample number: 1
------Sample number: 0
--------Sample number: 0
----------Sample number: 1
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
#ifndef THREAD_NUM
#    define THREAD_NUM   4
#endif
#ifndef BIN_SIZE
#    define BIN_SIZE   4
#endif
using namespace std;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( int rt, int val ) {
    if ( val <= 0) {
;        return;
    }
    if (RT.find(rt) == RT.end()) { 
        RT[rt] = val;
    } else {
        RT[rt] += val;
    }
    return;
}
void subBlkRT(int rt) {
    int msb = 0;
    int tmp_rt = rt;
    while(tmp_rt != 0) {
        tmp_rt = tmp_rt / 2;
        ++msb;
    }
    if (msb >= BIN_SIZE) {
        int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;
        for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {
            if (rt < b) {
                rtHistoCal(b - diff, 1);
                break;
            }
        }
    }
    else {
        rtHistoCal(pow(2, msb-1), 1);
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
/* A_addr (((((i + 1) * 256) * 256) + (j * 256)) + k) 0 */
int calAddrA_addr0( int t, int i, int j, int k) {
    int result = ((((((i + 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 1 */
int calAddrA_addr1( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i - 1) * 256) * 256) + (j * 256)) + k) 2 */
int calAddrA_addr2( int t, int i, int j, int k) {
    int result = ((((((i - 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + ((j + 1) * 256)) + k) 3 */
int calAddrA_addr3( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j + 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 4 */
int calAddrA_addr4( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + ((j - 1) * 256)) + k) 5 */
int calAddrA_addr5( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j - 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i * 256) * 256) + (j * 256)) + k) + 1) 6 */
int calAddrA_addr6( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) + 1)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 7 */
int calAddrA_addr7( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr (((((i * 256) * 256) + (j * 256)) + k) - 1) 8 */
int calAddrA_addr8( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) - 1)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 9 */
int calAddrA_addr9( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 0 */
int calAddrB_addr0( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i + 1) * 256) * 256) + (j * 256)) + k) 1 */
int calAddrB_addr1( int t, int i, int j, int k) {
    int result = ((((((i + 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 2 */
int calAddrB_addr2( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i - 1) * 256) * 256) + (j * 256)) + k) 3 */
int calAddrB_addr3( int t, int i, int j, int k) {
    int result = ((((((i - 1) * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + ((j + 1) * 256)) + k) 4 */
int calAddrB_addr4( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j + 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 5 */
int calAddrB_addr5( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + ((j - 1) * 256)) + k) 6 */
int calAddrB_addr6( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + ((j - 1) * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i * 256) * 256) + (j * 256)) + k) + 1) 7 */
int calAddrB_addr7( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) + 1)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 8 */
int calAddrB_addr8( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* B_addr (((((i * 256) * 256) + (j * 256)) + k) - 1) 9 */
int calAddrB_addr9( int t, int i, int j, int k) {
    int result = ((((((i * 256) * 256) + (j * 256)) + k) - 1)) * 8 / 64;
    return result;
}
/* B_addr ((((i * 256) * 256) + (j * 256)) + k) 10 */
int calAddrB_addr10( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
/* A_addr ((((i * 256) * 256) + (j * 256)) + k) 10 */
int calAddrA_addr10( int t, int i, int j, int k) {
    int result = (((((i * 256) * 256) + (j * 256)) + k)) * 8 / 64;
    return result;
}
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr1(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr2(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr3(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr4(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 4 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr5(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 5 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr6() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr6(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 6 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr7() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr7(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 7 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr8() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr8(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 8 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 8 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr9() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 8 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr9(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 9 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            if ( t == t_Start ) {
                iLB1 = i_Start;
            }
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB3 = k_Start;
                    }
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            int threadLB = 0;
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    int threadLB = 0;
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 5 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 6 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 7 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 8 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 9 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr0(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 0 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_B_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 0 */
#endif
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 8 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 9 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrB_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrB_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrB_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr10() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 1 + 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (255 - 1) + 1;
        if ( (255 - 1) == 0) goto SAMPLE;
        int k_Start = rand() % (255 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" + std::to_string(k_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 4  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDV: k  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 1) / THREAD_NUM;
        seperator = (10 - 1) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 1 + (i) * (B+1);
                BLIST[i][1] = 1 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 1 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 1 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
            }
        }
        thread_Start = 0;
        {
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (t_Start >= BLIST[i][0] && t_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int tLB0 = t_Start;
        int threadLB = 0;
        for ( int t = tLB0; t <= 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB1 = 1;
            int threadLB = 0;
            for ( int i = iLB1; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB2 = 1;
                int threadLB = 0;
                for ( int j = jLB2; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB3 = 1;
                    int threadLB = 0;
                    for ( int k = kLB3; k < 255; k++) {
                        if ( t == t_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr1( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr1(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr2( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr2(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 2 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr3( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr3(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 3 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr4( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr4(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 4 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr5( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr5(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 5 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr6( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr6(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 6 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr7( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr7(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 7 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr8( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr8(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 8 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr9( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr9(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 9 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 1;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 255; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 255; j++) {
                    /* Generating thread local iteration space mapping code */
                    {
                    int kLB6 = 1;
                    if ( t == t_Start && i == i_Start && j == j_Start ) {
                        kLB6 = k_Start;
                    }
                    for ( int k = kLB6; k < 255; k++) {
                        if ( t == t_Start  && i == i_Start  && j == j_Start && k == k_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = t - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(B_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << ", "<< k << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_t = BLIST[tid][0] + gap;
                            if (thread_t > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_i = i;
                            int thread_j = j;
                            int thread_k = k;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr10( thread_t, thread_i, thread_j, thread_k) == calAddrA_addr10(t_Start, i_Start, j_Start, k_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr10(thread_t, thread_i, thread_j, thread_k) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", "<< thread_k << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 10 refNumber[LoopRefTree->AA]: 10 */
#endif
                        threadLB = 0;
                    } // end of outer for loops
                    }
                    if (cntStart == true) {
                        threadLB = 0;
                    }
                } // end of outer for loops
                }
                if (cntStart == true) {
                    threadLB = 0;
                }
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
int main() {
    ref_A_addr3();
    ref_A_addr0();
    ref_A_addr5();
    ref_A_addr1();
    ref_A_addr2();
    ref_A_addr7();
    ref_A_addr8();
    ref_A_addr4();
    ref_B_addr1();
    ref_B_addr2();
    ref_B_addr3();
    ref_A_addr6();
    ref_B_addr4();
    ref_B_addr5();
    ref_B_addr6();
    ref_B_addr7();
    ref_B_addr8();
    ref_B_addr9();
    ref_A_addr9();
    ref_B_addr0();
    ref_B_addr10();
    ref_A_addr10();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: heat_3d */ 
