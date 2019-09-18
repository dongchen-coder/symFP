
 /* Start to analysis array index
Array index info: Total number of references: 20
y1.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgIn.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
y1.addr ((i * 1024) + j)
y2.addr ((i * 1024) + j)
imgOut.addr ((i * 1024) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %y1
double* %imgIn
double* %y2
double* %imgOut
double %alpha

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
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (1023, 0)
----Loop inc: (j + -1)
----Loop predicate: >=
------array access y2.addr ((i * 1024) + j)
------array access imgIn.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (0, 1024)
----Loop inc: (i + 1)
----Loop predicate: <
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y1.addr ((i * 1024) + j)
--j
--Loop Bound: (0, 1024)
--Loop inc: (j + 1)
--Loop predicate: <
----i
----Loop Bound: (1023, 0)
----Loop inc: (i + -1)
----Loop predicate: >=
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access y1.addr ((i * 1024) + j)
------array access y2.addr ((i * 1024) + j)
------array access imgOut.addr ((i * 1024) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
----Sample number: 10
------Sample number: 104
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
/* imgIn_addr ((i * 1024) + j) 0 */
int calAddrimgIn_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 0 */
int calAddry1_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgIn_addr ((i * 1024) + j) 1 */
int calAddrimgIn_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 1 */
int calAddry1_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 0 */
int calAddry2_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgIn_addr ((i * 1024) + j) 2 */
int calAddrimgIn_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 1 */
int calAddry2_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 2 */
int calAddry1_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 2 */
int calAddry2_addr2( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 0 */
int calAddrimgOut_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 1 */
int calAddrimgOut_addr1( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 3 */
int calAddry1_addr3( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 2 */
int calAddrimgOut_addr2( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 4 */
int calAddry1_addr4( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 3 */
int calAddry2_addr3( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 3 */
int calAddrimgOut_addr3( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 4 */
int calAddry2_addr4( int j, int i) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y1_addr ((i * 1024) + j) 5 */
int calAddry1_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* y2_addr ((i * 1024) + j) 5 */
int calAddry2_addr5( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* imgOut_addr ((i * 1024) + j) 4 */
int calAddrimgOut_addr4( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_y1_addr0() {
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry1_addr0( thread_i, thread_j) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry1_addr1( thread_i, thread_j) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
            int jLB3 = 1023;
            int threadLB = 0;
            for ( int j = jLB3; j >= 0; j--) {
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                    if ( calAddry1_addr2( thread_i, thread_j) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( thread_j, thread_i) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr4( thread_j, thread_i) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgIn_addr1() {
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
                    if ( calAddrimgIn_addr0( thread_i, thread_j) == calAddrimgIn_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgIn_addr1( thread_i, thread_j) == calAddrimgIn_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            int jLB3 = 1023;
            int threadLB = 0;
            for ( int j = jLB3; j >= 0; j--) {
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgIn_addr2( thread_i, thread_j) == calAddrimgIn_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y1_addr1() {
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry1_addr0( thread_i, thread_j) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry1_addr1( thread_i, thread_j) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
            int jLB3 = 1023;
            int threadLB = 0;
            for ( int j = jLB3; j >= 0; j--) {
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                    if ( calAddry1_addr2( thread_i, thread_j) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( thread_j, thread_i) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr4( thread_j, thread_i) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y2_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
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
                    if ( calAddry2_addr0( thread_i, thread_j) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr1( thread_i, thread_j) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr2( thread_i, thread_j) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr3( thread_j, thread_i) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( thread_j, thread_i) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr1(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y1_addr2() {
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
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
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
                    if ( calAddry1_addr2( thread_i, thread_j) == calAddry1_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( thread_j, thread_i) == calAddry1_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr4( thread_j, thread_i) == calAddry1_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgIn_addr0() {
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
                    if ( calAddrimgIn_addr0( thread_i, thread_j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgIn_addr1( thread_i, thread_j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            int jLB3 = 1023;
            int threadLB = 0;
            for ( int j = jLB3; j >= 0; j--) {
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgIn_addr2( thread_i, thread_j) == calAddrimgIn_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y2_addr2() {
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
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr2( thread_i, thread_j) == calAddry2_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr3( thread_j, thread_i) == calAddry2_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( thread_j, thread_i) == calAddry2_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgOut_addr0() {
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
        int iLB4 = i_Start;
        int threadLB = 0;
        for ( int i = iLB4; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            if ( i == i_Start ) {
                jLB5 = j_Start;
            }
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr0( thread_i, thread_j) == calAddrimgOut_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr1( thread_j, thread_i) == calAddrimgOut_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr1(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr2( thread_j, thread_i) == calAddrimgOut_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr2(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr3( thread_j, thread_i) == calAddrimgOut_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr4( thread_i, thread_j) == calAddrimgOut_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
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
void ref_y1_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB6 = j_Start;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( thread_j, thread_i) == calAddry1_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr4( thread_j, thread_i) == calAddry1_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgOut_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB6 = j_Start;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr1( thread_j, thread_i) == calAddrimgOut_addr2(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr1(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr2( thread_j, thread_i) == calAddrimgOut_addr2(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr2(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr3( thread_j, thread_i) == calAddrimgOut_addr2(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr4( thread_i, thread_j) == calAddrimgOut_addr2(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 4 */
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
void ref_y1_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB6 = j_Start;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr3( thread_j, thread_i) == calAddry1_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr4( thread_j, thread_i) == calAddry1_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y2_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
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
                    if ( calAddry2_addr0( thread_i, thread_j) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(imgIn_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr1( thread_i, thread_j) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr2( thread_i, thread_j) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr3( thread_j, thread_i) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( thread_j, thread_i) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 0 refNumber[LoopRefTree->AA]: 4 */
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr0(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgIn_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int i_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int j_Start = rand() % (1023 - 0 + 1) + 0;
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
            int jLB3 = 1023;
            if ( i == i_Start ) {
                jLB3 = j_Start;
            }
            for ( int j = jLB3; j >= 0; j--) {
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgIn_addr2( thread_i, thread_j) == calAddrimgIn_addr2(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgIn_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
            /* Generating thread local iteration space mapping code */
            {
            int jLB5 = 0;
            int threadLB = 0;
            for ( int j = jLB5; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
        int jLB6 = 0;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            int threadLB = 0;
            for ( int i = iLB7; i < 1024; i++) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = i - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgOut_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (1024 - 0) == 0) goto SAMPLE;
        int i_Start = rand() % (1024 - 0) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB6 = j_Start;
        int threadLB = 0;
        for ( int j = jLB6; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB7 = 0;
            if ( j == j_Start ) {
                iLB7 = i_Start;
            }
            for ( int i = iLB7; i < 1024; i++) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr1( thread_j, thread_i) == calAddrimgOut_addr1(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr1(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr2( thread_j, thread_i) == calAddrimgOut_addr1(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr2(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y1_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int jLB8 = 0;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            int threadLB = 0;
            for ( int i = iLB9; i >= 0; i--) {
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr3( thread_j, thread_i) == calAddrimgOut_addr1(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << ", " << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr4( thread_i, thread_j) == calAddrimgOut_addr1(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 1 refNumber[LoopRefTree->AA]: 4 */
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
void ref_y2_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB8 = j_Start;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr3( thread_j, thread_i) == calAddry2_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( thread_j, thread_i) == calAddry2_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgOut_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB8 = j_Start;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddrimgOut_addr3( thread_j, thread_i) == calAddrimgOut_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
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
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr4( thread_i, thread_j) == calAddrimgOut_addr3(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 3 refNumber[LoopRefTree->AA]: 4 */
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
void ref_y2_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 104;) {
SAMPLE:
        int j_Start = rand() % (1024 - 0) + 0;
        if ( (0 - 1023 + 1) == 0) goto SAMPLE;
        int i_Start = rand() % (1023 - 0 + 1) + 0;
        string idx_string = std::to_string(j_Start) + "_" + std::to_string(i_Start) + "_" ;
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
        /* Sampled IDV: j  */
        /* Sampled IDV: i  */
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
            if (j_Start >= BLIST[i][0] && j_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int jLB8 = j_Start;
        int threadLB = 0;
        for ( int j = jLB8; j < 1024; j++) {
            /* Generating thread local iteration space mapping code */
            {
            int iLB9 = 1023;
            if ( j == j_Start ) {
                iLB9 = i_Start;
            }
            for ( int i = iLB9; i >= 0; i--) {
                if ( j == j_Start && i == i_Start && !cntStart ) {
                    threadLB = thread_Start;
                }
                int gap = j - BLIST[thread_Start][0];
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr3( thread_j, thread_i) == calAddry2_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr3(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
#endif
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << j << ", "<< i << "), " << cnt << ") " << endl;
#endif
                    }
                } // end of interleaving loop
                /* Generating thread local iteration space mapping code */
                for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                    int thread_j = BLIST[tid][0] + gap;
                    if (thread_j > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_i = i;
                if (cntStart == true) {
                    cnt++;
                    if ( calAddry2_addr4( thread_j, thread_i) == calAddry2_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr4(thread_j, thread_i) << ", " << "(" << thread_j << ", "<< thread_i << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
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
        int iLB10 = 0;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            int threadLB = 0;
            for ( int j = jLB11; j < 1024; j++) {
                int gap = j - BLIST[thread_Start][0];
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr4(j_Start, i_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << ", " << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y1_addr5() {
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
        int iLB10 = i_Start;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
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
                    if ( calAddry1_addr5( thread_i, thread_j) == calAddry1_addr5(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry1_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_y2_addr5() {
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
        int iLB10 = i_Start;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddry2_addr5( thread_i, thread_j) == calAddry2_addr5(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddry2_addr5(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_i = BLIST[tid][0] + gap;
                    if (thread_i > BLIST[tid][1]) {
                        continue;
                    }
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
#ifdef DEBUG
                        cout << "(imgOut_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_imgOut_addr4() {
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
        int iLB10 = i_Start;
        int threadLB = 0;
        for ( int i = iLB10; i < 1024; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB11 = 0;
            if ( i == i_Start ) {
                jLB11 = j_Start;
            }
            for ( int j = jLB11; j < 1024; j++) {
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
                        cout << "(y1_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(y2_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                    if ( calAddrimgOut_addr4( thread_i, thread_j) == calAddrimgOut_addr4(i_Start, j_Start)) {
                            subBlkRT(cnt);
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrimgOut_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 4 */
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
    ref_y1_addr0();
    ref_imgIn_addr1();
    ref_y1_addr1();
    ref_y2_addr1();
    ref_y1_addr2();
    ref_imgIn_addr0();
    ref_y2_addr2();
    ref_imgOut_addr0();
    ref_y1_addr3();
    ref_imgOut_addr2();
    ref_y1_addr4();
    ref_y2_addr0();
    ref_imgIn_addr2();
    ref_imgOut_addr1();
    ref_y2_addr3();
    ref_imgOut_addr3();
    ref_y2_addr4();
    ref_y1_addr5();
    ref_y2_addr5();
    ref_imgOut_addr4();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: deriche */ 
