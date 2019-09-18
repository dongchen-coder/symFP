
 /* Start to analysis array index
Array index info: Total number of references: 16
_fict_.addr t
ey.addr ((i * 1024) + j)
ey.addr (0 + j)
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)
hz.addr (((i * 1024) + j) - 1)
hz.addr ((i * 1024) + j)
hz.addr (((i - 1) * 1024) + j)
hz.addr ((i * 1024) + j)
ex.addr (((i * 1024) + j) + 1)
ex.addr ((i * 1024) + j)
ex.addr ((i * 1024) + j)
ex.addr ((i * 1024) + j)
ey.addr (((i + 1) * 1024) + j)
ey.addr ((i * 1024) + j)
hz.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %_fict_
double* %ey
double* %ex
double* %hz

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--t
--Loop Bound: (0, 10)
--Loop inc: (t + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access _fict_.addr t
------array access ey.addr (0 + j)
----i
----Loop Bound: (1, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 1024)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i - 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (1, 1024)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access ex.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)
--------array access hz.addr (((i * 1024) + j) - 1)
--------array access ex.addr ((i * 1024) + j)
----i
----Loop Bound: (0, 1023)
----Loop inc: (i + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 1023)
------Loop inc: (j + 1)
------Loop predicate: <
--------array access hz.addr ((i * 1024) + j)
--------array access ex.addr (((i * 1024) + j) + 1)
--------array access ex.addr ((i * 1024) + j)
--------array access ey.addr (((i + 1) * 1024) + j)
--------array access ey.addr ((i * 1024) + j)
--------array access hz.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 0
------Sample number: 1
------Sample number: 1
--------Sample number: 10
------Sample number: 1
--------Sample number: 10
------Sample number: 1
--------Sample number: 10
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
/* _fict__addr t 0 */
int calAddr_fict__addr0( int t, int j) {
    int result = (t) * 8 / 64;
    return result;
}
/* ey_addr (0 + j) 0 */
int calAddrey_addr0( int t, int j) {
    int result = ((0 + j)) * 8 / 64;
    return result;
}
/* ey_addr ((i * 1024) + j) 1 */
int calAddrey_addr1( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr ((i * 1024) + j) 0 */
int calAddrhz_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr (((i - 1) * 1024) + j) 1 */
int calAddrhz_addr1( int t, int i, int j) {
    int result = ((((i - 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* ey_addr ((i * 1024) + j) 2 */
int calAddrey_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* ex_addr ((i * 1024) + j) 0 */
int calAddrex_addr0( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr ((i * 1024) + j) 2 */
int calAddrhz_addr2( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr (((i * 1024) + j) - 1) 3 */
int calAddrhz_addr3( int t, int i, int j) {
    int result = ((((i * 1024) + j) - 1)) * 8 / 64;
    return result;
}
/* ex_addr ((i * 1024) + j) 1 */
int calAddrex_addr1( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr ((i * 1024) + j) 4 */
int calAddrhz_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* ex_addr (((i * 1024) + j) + 1) 2 */
int calAddrex_addr2( int t, int i, int j) {
    int result = ((((i * 1024) + j) + 1)) * 8 / 64;
    return result;
}
/* ex_addr ((i * 1024) + j) 3 */
int calAddrex_addr3( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* ey_addr (((i + 1) * 1024) + j) 3 */
int calAddrey_addr3( int t, int i, int j) {
    int result = ((((i + 1) * 1024) + j)) * 8 / 64;
    return result;
}
/* ey_addr ((i * 1024) + j) 4 */
int calAddrey_addr4( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* hz_addr ((i * 1024) + j) 5 */
int calAddrhz_addr5( int t, int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref__fict__addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 2  */
        /* Sampled IDV: t  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( t == t_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( t == t_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = t - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddr_fict__addr0( thread_t, thread_j) == calAddr_fict__addr0(t_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddr_fict__addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ey_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( thread_t, thread_j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrey_addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( thread_t, thread_i, thread_j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( thread_t, thread_i, thread_j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 2 refNumber[LoopRefTree->AA]: 2 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( thread_t, thread_i, thread_j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( thread_t, thread_i, thread_j) == calAddrey_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ey_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 1;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 2  */
        /* Sampled IDV: t  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( t == t_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( t == t_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = t - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( thread_t, thread_j) == calAddrey_addr0(t_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrey_addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( thread_t, thread_i, thread_j) == calAddrey_addr0(t_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( thread_t, thread_i, thread_j) == calAddrey_addr0(t_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 0 refNumber[LoopRefTree->AA]: 2 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( thread_t, thread_i, thread_j) == calAddrey_addr0(t_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( thread_t, thread_i, thread_j) == calAddrey_addr0(t_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ey_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( thread_t, thread_j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrey_addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( thread_t, thread_i, thread_j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( thread_t, thread_i, thread_j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 1 refNumber[LoopRefTree->AA]: 2 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( thread_t, thread_i, thread_j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( thread_t, thread_i, thread_j) == calAddrey_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_hz_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 2 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
void ref_hz_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 3 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
void ref_hz_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 0 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
void ref_hz_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 1) + 1;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            if ( t == t_Start ) {
                iLB2 = i_Start;
            }
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB3 = j_Start;
                }
                for ( int j = jLB3; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 1 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
void ref_hz_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 4 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
void ref_ex_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( thread_t, thread_i, thread_j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( thread_t, thread_i, thread_j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( thread_t, thread_i, thread_j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( thread_t, thread_i, thread_j) == calAddrex_addr2(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ex_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( thread_t, thread_i, thread_j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( thread_t, thread_i, thread_j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( thread_t, thread_i, thread_j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( thread_t, thread_i, thread_j) == calAddrex_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ex_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( thread_t, thread_i, thread_j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( thread_t, thread_i, thread_j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 0 refNumber[LoopRefTree->AA]: 1 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( thread_t, thread_i, thread_j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( thread_t, thread_i, thread_j) == calAddrex_addr0(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ex_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 1) + 1;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            if ( t == t_Start ) {
                iLB4 = i_Start;
            }
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                if ( t == t_Start && i == i_Start ) {
                    jLB5 = j_Start;
                }
                for ( int j = jLB5; j < 1024; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr0( thread_t, thread_i, thread_j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr1( thread_t, thread_i, thread_j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 1 refNumber[LoopRefTree->AA]: 1 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            int threadLB = 0;
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                int threadLB = 0;
                for ( int j = jLB7; j < 1023; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr2( thread_t, thread_i, thread_j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrex_addr3( thread_t, thread_i, thread_j) == calAddrex_addr1(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrex_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ey_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( thread_t, thread_j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrey_addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( thread_t, thread_i, thread_j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( thread_t, thread_i, thread_j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 3 refNumber[LoopRefTree->AA]: 2 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( thread_t, thread_i, thread_j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( thread_t, thread_i, thread_j) == calAddrey_addr3(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_ey_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrey_addr0( thread_t, thread_j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrey_addr0(thread_t, thread_j) << ", " << "(" << thread_t << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr1( thread_t, thread_i, thread_j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr2( thread_t, thread_i, thread_j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 4 refNumber[LoopRefTree->AA]: 2 */
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr3( thread_t, thread_i, thread_j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrey_addr4( thread_t, thread_i, thread_j) == calAddrey_addr4(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrey_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(hz_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_hz_addr5() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 10;) {
SAMPLE:
        int t_Start = rand() % (10 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0) + 0;
        if ( (1023 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0) + 0;
        string idx_string = std::to_string(t_Start) + "_" + std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDVs 3  */
        /* Sampled IDV: t  */
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (10 - 0) / THREAD_NUM;
        seperator = (10 - 0) - THREAD_NUM * B;
        for (int i = 0; i < THREAD_NUM; ++i) {
            if (i < seperator) {
                BLIST[i][0] = 0 + (i) * (B+1);
                BLIST[i][1] = 0 + (i+1) * (B+1) - 1;
            }
            else {
                BLIST[i][0] = 0 + seperator * (B+1) + (i - seperator) * B;
                BLIST[i][1] = 0 + seperator * (B+1) + (i - seperator + 1)  * B - 1;
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
        for ( int t = tLB0; t < 10; t++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            int threadLB = 0;
            for ( int j = jLB1; j < 1024; j++) {
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(_fict__addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_t = BLIST[tid][0] + gap;
                    if (thread_t > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(ey_addr " <<  ", " << "(" << t << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
            } // end of outer for loops
            }
            if (cntStart == true) {
                threadLB = 0;
            }
            /* Generating thread local iteration space mapping code */
            {
            int iLB2 = 1;
            int threadLB = 0;
            for ( int i = iLB2; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB3 = 0;
                int threadLB = 0;
                for ( int j = jLB3; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr0( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr0(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr1( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr1(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB4 = 0;
            int threadLB = 0;
            for ( int i = iLB4; i < 1024; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB5 = 1;
                int threadLB = 0;
                for ( int j = jLB5; j < 1024; j++) {
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr2( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr2(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr3( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr3(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int iLB6 = 0;
            if ( t == t_Start ) {
                iLB6 = i_Start;
            }
            for ( int i = iLB6; i < 1023; i++) {
                /* Generating thread local iteration space mapping code */
                {
                int jLB7 = 0;
                if ( t == t_Start && i == i_Start ) {
                    jLB7 = j_Start;
                }
                for ( int j = jLB7; j < 1023; j++) {
                    if ( t == t_Start  && i == i_Start && j == j_Start && !cntStart ) {
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr4( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr4(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ex_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(ey_addr " <<  ", " << "(" << t << ", "<< i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrhz_addr5( thread_t, thread_i, thread_j) == calAddrhz_addr5(t_Start, i_Start, j_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrhz_addr5(thread_t, thread_i, thread_j) << ", " << "(" << thread_t << ", "<< thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 5 refNumber[LoopRefTree->AA]: 5 */
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
EndSample:
        s++;
        }
}
int main() {
    ref__fict__addr0();
    ref_ey_addr2();
    ref_ey_addr0();
    ref_ey_addr1();
    ref_hz_addr2();
    ref_hz_addr3();
    ref_hz_addr0();
    ref_hz_addr1();
    ref_hz_addr4();
    ref_ex_addr2();
    ref_ex_addr3();
    ref_ex_addr0();
    ref_ex_addr1();
    ref_ey_addr3();
    ref_ey_addr4();
    ref_hz_addr5();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: fdtd_2d */ 
