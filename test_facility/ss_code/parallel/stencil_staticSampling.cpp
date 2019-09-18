
 /* Start to analysis array index
Array index info: Total number of references: 6
a.addr ((6 * i) + j)
a.addr (((6 * i) + j) + 1)
b.addr ((6 * i) + j)
a.addr (((6 * i) + j) - 1)
a.addr ((6 * (i - 1)) + j)
a.addr ((6 * (i + 1)) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %a
double* %b

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (1, 5)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (1, 5)
----Loop inc: (j + 1)
----Loop predicate: <
------array access a.addr ((6 * i) + j)
------array access a.addr (((6 * i) + j) + 1)
------array access a.addr (((6 * i) + j) - 1)
------array access a.addr ((6 * (i - 1)) + j)
------array access a.addr ((6 * (i + 1)) + j)
------array access b.addr ((6 * i) + j)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 4
------Sample number: 16
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
/* a_addr ((6 * i) + j) 0 */
int calAddra_addr0( int i, int j) {
    int result = (((6 * i) + j)) * 8 / 64;
    return result;
}
/* a_addr (((6 * i) + j) + 1) 1 */
int calAddra_addr1( int i, int j) {
    int result = ((((6 * i) + j) + 1)) * 8 / 64;
    return result;
}
/* a_addr (((6 * i) + j) - 1) 2 */
int calAddra_addr2( int i, int j) {
    int result = ((((6 * i) + j) - 1)) * 8 / 64;
    return result;
}
/* a_addr ((6 * (i - 1)) + j) 3 */
int calAddra_addr3( int i, int j) {
    int result = (((6 * (i - 1)) + j)) * 8 / 64;
    return result;
}
/* a_addr ((6 * (i + 1)) + j) 4 */
int calAddra_addr4( int i, int j) {
    int result = (((6 * (i + 1)) + j)) * 8 / 64;
    return result;
}
/* b_addr ((6 * i) + j) 0 */
int calAddrb_addr0( int i, int j) {
    int result = (((6 * i) + j)) * 8 / 64;
    return result;
}
void ref_a_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        if ( calAddra_addr0( thread_i, thread_j) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr1( thread_i, thread_j) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr2( thread_i, thread_j) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr3( thread_i, thread_j) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr4( thread_i, thread_j) == calAddra_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(b_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_a_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        if ( calAddra_addr0( thread_i, thread_j) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr1( thread_i, thread_j) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr2( thread_i, thread_j) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr3( thread_i, thread_j) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr4( thread_i, thread_j) == calAddra_addr1(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(b_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_b_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        cout << "(a_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(a_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(a_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(a_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        cout << "(a_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
                        if ( calAddrb_addr0( thread_i, thread_j) == calAddrb_addr0(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddrb_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
EndSample:
        s++;
        }
}
void ref_a_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        if ( calAddra_addr0( thread_i, thread_j) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
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
                        if ( calAddra_addr1( thread_i, thread_j) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
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
                        if ( calAddra_addr2( thread_i, thread_j) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr3( thread_i, thread_j) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr4( thread_i, thread_j) == calAddra_addr2(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(b_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_a_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        if ( calAddra_addr0( thread_i, thread_j) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 3 refNumber[LoopRefTree->AA]: 0 */
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
                        if ( calAddra_addr1( thread_i, thread_j) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
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
                        if ( calAddra_addr2( thread_i, thread_j) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                    int thread_j = j;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddra_addr3( thread_i, thread_j) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        if ( calAddra_addr4( thread_i, thread_j) == calAddra_addr3(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(b_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
void ref_a_addr4() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int i_Start = rand() % (5 - 1) + 1;
        if ( (5 - 1) == 0) goto SAMPLE;
        int j_Start = rand() % (5 - 1) + 1;
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
        B = (5 - 1) / THREAD_NUM;
        seperator = (5 - 1) - THREAD_NUM * B;
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
            if (i_Start >= BLIST[i][0] && i_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int iLB0 = i_Start;
        int threadLB = 0;
        for ( int i = iLB0; i < 5; i++) {
            /* Generating thread local iteration space mapping code */
            {
            int jLB1 = 1;
            if ( i == i_Start ) {
                jLB1 = j_Start;
            }
            for ( int j = jLB1; j < 5; j++) {
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
                        if ( calAddra_addr0( thread_i, thread_j) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr0(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 0 */
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
                        if ( calAddra_addr1( thread_i, thread_j) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr1(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 1 */
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
                        if ( calAddra_addr2( thread_i, thread_j) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr2(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 2 */
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
                        if ( calAddra_addr3( thread_i, thread_j) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr3(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
                    }
#endif
                } // end of interleaving loop
#ifdef DEBUG
                cout << endl;
                /* useID: 4 refNumber[LoopRefTree->AA]: 3 */
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
                        if ( calAddra_addr4( thread_i, thread_j) == calAddra_addr4(i_Start, j_Start)) {
#ifdef DEBUG
                            rtHistoCal(cnt, 1);
#else
                            subBlkRT(cnt);
#endif
                        goto EndSample;
                    }
                }
                    cntStart = true;
#ifdef DEBUG
                    if (cntStart == true) {
                        cout << "(" << calAddra_addr4(thread_i, thread_j) << ", " << "(" << thread_i << ", "<< thread_j << "), " << cnt << ") " << endl;
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
                        cout << "(b_addr " <<  ", " << "(" << i << ", "<< j << "), " << cnt << ") " << endl;
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
int main() {
    ref_a_addr0();
    ref_a_addr1();
    ref_b_addr0();
    ref_a_addr2();
    ref_a_addr3();
    ref_a_addr4();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: stencil */ 
