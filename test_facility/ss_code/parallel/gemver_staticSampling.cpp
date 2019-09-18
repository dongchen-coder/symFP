
 /* Start to analysis array index
Array index info: Total number of references: 17
A.addr ((i * 1024) + j)
A.addr ((i * 1024) + j)
u1.addr i
y.addr j
x.addr i
v1.addr j
u2.addr i
v2.addr j
x.addr i
w.addr i
x.addr i
A.addr ((j * 1024) + i)
x.addr i
z.addr i
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
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access A.addr ((i * 1024) + j)
------array access u1.addr i
------array access v1.addr j
------array access u2.addr i
------array access v2.addr j
------array access A.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access x.addr i
------array access A.addr ((j * 1024) + i)
------array access y.addr j
------array access x.addr i
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----array access x.addr i
----array access z.addr i
----array access x.addr i
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access w.addr i
------array access A.addr ((i * 1024) + j)
------array access x.addr j
------array access w.addr i

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
----Sample number: 10
------Sample number: 104
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
/* A_addr ((i * 1024) + j) 0 */
int calAddrA_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* u1_addr i 0 */
int calAddru1_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* v1_addr j 0 */
int calAddrv1_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* u2_addr i 0 */
int calAddru2_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* v2_addr j 0 */
int calAddrv2_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 1 */
int calAddrA_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* x_addr i 0 */
int calAddrx_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* A_addr ((j * 1024) + i) 2 */
int calAddrA_addr2( int i, int j) {
    int result = (((j * 1024) + i)) * 8 / 64;
    return result;
}
/* y_addr j 0 */
int calAddry_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* x_addr i 1 */
int calAddrx_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* x_addr i 2 */
int calAddrx_addr2( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* z_addr i 0 */
int calAddrz_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* x_addr i 3 */
int calAddrx_addr3( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* w_addr i 0 */
int calAddrw_addr0( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + j) 3 */
int calAddrA_addr3( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* x_addr j 4 */
int calAddrx_addr4( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* w_addr i 1 */
int calAddrw_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( thread_i, thread_j) == calAddrA_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( thread_i, thread_j) == calAddrA_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( thread_i, thread_j) == calAddrA_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( thread_i, thread_j) == calAddrA_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr0( thread_i, thread_j) == calAddrA_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr1( thread_i, thread_j) == calAddrA_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( thread_i, thread_j) == calAddrA_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( thread_i, thread_j) == calAddrA_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_u1_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru1_addr0( thread_i, thread_j) == calAddru1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddru1_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_y_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB2 = i_Start;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry_addr0( thread_i, thread_j) == calAddry_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB2 = i_Start;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( thread_i, thread_j) == calAddrx_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( thread_i, thread_j) == calAddrx_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( thread_i) == calAddrx_addr1(i_Start, j_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr2(thread_i) << ", " << "(" << thread_i << ", " << "), " << cnt << ") " << endl;
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
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( thread_i) == calAddrx_addr1(i_Start, j_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr3(thread_i) << ", " << "(" << thread_i << ", " << "), " << cnt << ") " << endl;
                }
#endif
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
            /* useID: 1 refNumber[LoopRefTree->AA]: 3 */
#endif
            threadLB = 0;
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( thread_i, thread_j) == calAddrx_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_v1_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv1_addr0( thread_i, thread_j) == calAddrv1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrv1_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_u2_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddru2_addr0( thread_i, thread_j) == calAddru2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddru2_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_v2_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 0;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(v1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(u2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrv2_addr0( thread_i, thread_j) == calAddrv2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrv2_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB2 = 0;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            int threadLB = 0;
            for ( int j = jLB3; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr3() {
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

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 1  */
        /* Sampled IDV: i  */
        /* Sampled IDVs 1  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if ( i == i_Start && !cntStart ) {
                threadLB = thread_Start;
            }
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( thread_i) == calAddrx_addr3(i_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr2(thread_i) << ", " << "(" << thread_i << "), " << cnt << ") " << endl;
                }
#endif
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
            /* useID: 3 refNumber[LoopRefTree->AA]: 2 */
#endif
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( thread_i) == calAddrx_addr3(i_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
                cntStart = true;
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr3(thread_i) << ", " << "(" << thread_i << "), " << cnt << ") " << endl;
                }
#endif
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
            /* useID: 3 refNumber[LoopRefTree->AA]: 3 */
#endif
            threadLB = 0;
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( thread_i, thread_j) == calAddrx_addr3(i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_w_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB5 = i_Start;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( thread_i, thread_j) == calAddrw_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrw_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( thread_i, thread_j) == calAddrw_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrw_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB2 = i_Start;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr0( thread_i, thread_j) == calAddrx_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr1( thread_i, thread_j) == calAddrx_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( thread_i) == calAddrx_addr0(i_Start, j_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr2(thread_i) << ", " << "(" << thread_i << ", " << "), " << cnt << ") " << endl;
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
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( thread_i) == calAddrx_addr0(i_Start, j_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr3(thread_i) << ", " << "(" << thread_i << ", " << "), " << cnt << ") " << endl;
                }
#endif
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
            /* useID: 0 refNumber[LoopRefTree->AA]: 3 */
#endif
            threadLB = 0;
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( thread_i, thread_j) == calAddrx_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_A_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB2 = i_Start;
        int threadLB = 0;
        for ( int i = iLB2; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB3 = 0;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr2( thread_i, thread_j) == calAddrA_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB4 = 0;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( thread_i, thread_j) == calAddrA_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr2() {
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

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 1  */
        /* Sampled IDV: i  */
        /* Sampled IDVs 1  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if ( i == i_Start && !cntStart ) {
                threadLB = thread_Start;
            }
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr2( thread_i) == calAddrx_addr2(i_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
                cntStart = true;
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr2(thread_i) << ", " << "(" << thread_i << "), " << cnt << ") " << endl;
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
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(z_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrx_addr3( thread_i) == calAddrx_addr2(i_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrx_addr3(thread_i) << ", " << "(" << thread_i << "), " << cnt << ") " << endl;
                }
#endif
            } // end of interleaving loop
#ifdef DEBUG
            cout << endl;
            /* useID: 2 refNumber[LoopRefTree->AA]: 3 */
#endif
            threadLB = 0;
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( thread_i, thread_j) == calAddrx_addr2(i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_z_addr0() {
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

        /* Variable used to compute thread-local iteration space */
        int B = 0;
        auto BLIST = new int[THREAD_NUM][2];
        int seperator = 0;
        int thread_Start = 0;
        /* Generating reuse search code */
        /* Sampled IDVs 1  */
        /* Sampled IDV: i  */
        /* Sampled IDVs 1  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            if ( i == i_Start && !cntStart ) {
                threadLB = thread_Start;
            }
            int gap = i - BLIST[thread_Start][0];
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
            /* Generating thread local iteration space mapping code */
            for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
            if (cntStart == true) {
                cnt++;
                if ( calAddrz_addr0( thread_i) == calAddrz_addr0(i_Start)) {
                        subBlkRT(cnt);
                    goto EndSample;
                }
            }
                cntStart = true;
#ifdef DEBUG
                if (cntStart == true) {
                    cout << "(" << calAddrz_addr0(thread_i) << ", " << "(" << thread_i << "), " << cnt << ") " << endl;
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
                int thread_i = BLIST[tid][0] + gap;
                if (thread_i > BLIST[tid][1]) {
                    continue;
                }
                if (cntStart == true) {
                    cnt++;
#ifdef DEBUG
                    cout << "(x_addr " <<  ", " << "(" << i << "), " << cnt << ") " << endl;
#endif
                }
            } // end of interleaving loop
        } // end of outer for loops
        }
        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
        int iLB5 = 0;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            int threadLB = 0;
            for ( int j = jLB6; j < 1024; j++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_A_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB5 = i_Start;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrA_addr3( thread_i, thread_j) == calAddrA_addr3(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrA_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_x_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB5 = i_Start;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrx_addr4( thread_i, thread_j) == calAddrx_addr4(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrx_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(w_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_w_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(i_Start) + "_" + std::to_string(j_Start) + "_" ;
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
        /* Sampled IDV: i  */
        /* Sampled IDV: j  */
        /* Sampled IDVs 2  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (1024 - 0) / THREAD_NUM;
        seperator = (1024 - 0) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB5 = i_Start;
        int threadLB = 0;
        for ( int i = iLB5; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB6 = 0;
            if ( i == i_Start ) {
                jLB6 = j_Start;
            }
            for ( int j = jLB6; j < 1024; j++) {
                if ( i == i_Start && j == j_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr0( thread_i, thread_j) == calAddrw_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrw_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(A_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(x_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrw_addr1( thread_i, thread_j) == calAddrw_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrw_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
int main() {
    ref_A_addr1();
    ref_A_addr0();
    ref_u1_addr0();
    ref_y_addr0();
    ref_x_addr1();
    ref_v1_addr0();
    ref_u2_addr0();
    ref_v2_addr0();
    ref_x_addr3();
    ref_w_addr0();
    ref_x_addr0();
    ref_A_addr2();
    ref_x_addr2();
    ref_z_addr0();
    ref_A_addr3();
    ref_x_addr4();
    ref_w_addr1();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: gemver */ 
