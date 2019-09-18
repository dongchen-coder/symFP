
 /* Start to analysis array index
Array index info: Total number of references: 7
sum.addr p
sum.addr p
A.addr ((((r * 256) * 256) + (q * 256)) + p)
A.addr ((((r * 256) * 256) + (q * 256)) + s)
C4.addr ((s * 256) + p)
sum.addr p
sum.addr p

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %sum
double* %A
double* %C4

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--r
--Loop Bound: (0, 256)
--Loop inc: (r + 1)
--Loop predicate: <
----q
----Loop Bound: (0, 256)
----Loop inc: (q + 1)
----Loop predicate: <
------p
------Loop Bound: (0, 256)
------Loop inc: (p + 1)
------Loop predicate: <
--------array access sum.addr p
--------s
--------Loop Bound: (0, 256)
--------Loop inc: (s + 1)
--------Loop predicate: <
----------array access A.addr ((((r * 256) * 256) + (q * 256)) + s)
----------array access C4.addr ((s * 256) + p)
----------array access sum.addr p
----------array access sum.addr p
------p
------Loop Bound: (0, 256)
------Loop inc: (p + 1)
------Loop predicate: <
--------array access sum.addr p
--------array access A.addr ((((r * 256) * 256) + (q * 256)) + p)

Finish analysis loops */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 2
------Sample number: 6
--------Sample number: 16
----------Sample number: 42
--------Sample number: 16
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
/* sum_addr p 0 */
int calAddrsum_addr0( int r, int q, int p) {
    int result = (p) * 8 / 64;
    return result;
}
/* A_addr ((((r * 256) * 256) + (q * 256)) + s) 0 */
int calAddrA_addr0( int r, int q, int p, int s) {
    int result = (((((r * 256) * 256) + (q * 256)) + s)) * 8 / 64;
    return result;
}
/* C4_addr ((s * 256) + p) 0 */
int calAddrC4_addr0( int r, int q, int p, int s) {
    int result = (((s * 256) + p)) * 8 / 64;
    return result;
}
/* sum_addr p 1 */
int calAddrsum_addr1( int r, int q, int p, int s) {
    int result = (p) * 8 / 64;
    return result;
}
/* sum_addr p 2 */
int calAddrsum_addr2( int r, int q, int p, int s) {
    int result = (p) * 8 / 64;
    return result;
}
/* sum_addr p 3 */
int calAddrsum_addr3( int r, int q, int p) {
    int result = (p) * 8 / 64;
    return result;
}
/* A_addr ((((r * 256) * 256) + (q * 256)) + p) 1 */
int calAddrA_addr1( int r, int q, int p) {
    int result = (((((r * 256) * 256) + (q * 256)) + p)) * 8 / 64;
    return result;
}
void ref_sum_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( thread_r, thread_q, thread_p) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr0(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 0 refNumber[LoopRefTree->AA]: 0 */
#endif
                    threadLB = 0;
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    int threadLB = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && p == p_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr1(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p<< thread_s << ", " << "), " << cnt << ") " << endl;
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
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr2(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p<< thread_s << ", " << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                int threadLB = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( thread_r, thread_q, thread_p) == calAddrsum_addr0(r_Start, q_Start, p_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr3(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
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
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
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
void ref_sum_addr3() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                int threadLB = 0;
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( thread_r, thread_q, thread_p) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr0(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 3 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    int threadLB = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr1(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 3 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr2(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << ", " << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB4 = p_Start;
                }
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( thread_r, thread_q, thread_p) == calAddrsum_addr3(r_Start, q_Start, p_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr3(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << "), " << cnt << ") " << endl;
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
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
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
void ref_A_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 16;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDVs 3  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                int threadLB = 0;
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    int threadLB = 0;
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && !cntStart) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_r, thread_q, thread_p, thread_s) == calAddrA_addr1(r_Start, q_Start, p_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << ", " << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB4 = p_Start;
                }
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( thread_r, thread_q, thread_p) == calAddrA_addr1(r_Start, q_Start, p_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
                        cntStart = true;
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrA_addr1(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << "), " << cnt << ") " << endl;
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
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDV: s  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && p == p_Start && s == s_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrA_addr0( thread_r, thread_q, thread_p, thread_s) == calAddrA_addr0(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrA_addr0(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
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
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                int threadLB = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr1( thread_r, thread_q, thread_p) == calAddrA_addr0(r_Start, q_Start, p_Start, s_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrA_addr1(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
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
        } // end of outer for loops
        }
EndSample:
        s++;
        }
}
void ref_C4_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDV: s  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && p == p_Start && s == s_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrC4_addr0( thread_r, thread_q, thread_p, thread_s) == calAddrC4_addr0(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrC4_addr0(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
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
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                int threadLB = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(sum_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
#endif
                        }
                    } // end of interleaving loop
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
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
void ref_sum_addr1() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDV: s  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( thread_r, thread_q, thread_p) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr0(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 1 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && p == p_Start && s == s_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr1(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
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
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr2(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                int threadLB = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( thread_r, thread_q, thread_p) == calAddrsum_addr1(r_Start, q_Start, p_Start, s_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr3(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
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
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
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
void ref_sum_addr2() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 42;) {
SAMPLE:
        int r_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int q_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int p_Start = rand() % (256 - 0) + 0;
        if ( (256 - 0) == 0) goto SAMPLE;
        int s_Start = rand() % (256 - 0) + 0;
        string idx_string = std::to_string(r_Start) + "_" + std::to_string(q_Start) + "_" + std::to_string(p_Start) + "_" + std::to_string(s_Start) + "_" ;
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
        /* Sampled IDV: r  */
        /* Sampled IDV: q  */
        /* Sampled IDV: p  */
        /* Sampled IDV: s  */
        /* Sampled IDVs 4  */

        /* Generating thread local iteration space mapping code */
#ifdef DEBUG
        cout << "Count: " << cnt << endl;
#endif
        B = (256 - 0) / THREAD_NUM;
        seperator = (256 - 0) - THREAD_NUM * B;
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
            if (r_Start >= BLIST[i][0] && r_Start <= BLIST[i][1] ) {
                thread_Start = i;
                break;
            }
        }
        int rLB0 = r_Start;
        int threadLB = 0;
        for ( int r = rLB0; r < 256; r++) {
            /* Generating thread local iteration space mapping code */
            {
            int qLB1 = 0;
            if ( r == r_Start ) {
                qLB1 = q_Start;
            }
            for ( int q = qLB1; q < 256; q++) {
                /* Generating thread local iteration space mapping code */
                {
                int pLB2 = 0;
                if ( r == r_Start && q == q_Start ) {
                    pLB2 = p_Start;
                }
                for ( int p = pLB2; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start && p == p_Start && !cntStart ) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr0( thread_r, thread_q, thread_p) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr0(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
                        }
#endif
                    } // end of interleaving loop
#ifdef DEBUG
                    cout << endl;
                    /* useID: 2 refNumber[LoopRefTree->AA]: 0 */
#endif
                    /* Generating thread local iteration space mapping code */
                    {
                    int sLB3 = 0;
                    if ( r == r_Start && q == q_Start && p == p_Start ) {
                        sLB3 = s_Start;
                    }
                    for ( int s = sLB3; s < 256; s++) {
                        if ( r == r_Start  && q == q_Start  && p == p_Start && s == s_Start && !cntStart ) {
                            threadLB = thread_Start;
                        }
                        int gap = r - BLIST[thread_Start][0];
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                            if (cntStart == true) {
                                cnt++;
#ifdef DEBUG
                                cout << "(C4_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << ", "<< s << "), " << cnt << ") " << endl;
#endif
                            }
                        } // end of interleaving loop
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr1( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr1(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
                            }
#endif
                        } // end of interleaving loop
#ifdef DEBUG
                        cout << endl;
                        /* useID: 2 refNumber[LoopRefTree->AA]: 1 */
#endif
                        /* Generating thread local iteration space mapping code */
                        for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                            int thread_r = BLIST[tid][0] + gap;
                            if (thread_r > BLIST[tid][1]) {
                                continue;
                            }
                            int thread_q = q;
                            int thread_p = p;
                            int thread_s = s;
                        if (cntStart == true) {
                            cnt++;
                            if ( calAddrsum_addr2( thread_r, thread_q, thread_p, thread_s) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                    subBlkRT(cnt);
                                goto EndSample;
                            }
                        }
                            cntStart = true;
#ifdef DEBUG
                            if (cntStart == true) {
                                cout << "(" << calAddrsum_addr2(thread_r, thread_q, thread_p, thread_s) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", "<< thread_s << "), " << cnt << ") " << endl;
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
                int pLB4 = 0;
                int threadLB = 0;
                for ( int p = pLB4; p < 256; p++) {
                    if ( r == r_Start  && q == q_Start  && !cntStart) {
                        threadLB = thread_Start;
                    }
                    int gap = r - BLIST[thread_Start][0];
                    /* Generating thread local iteration space mapping code */
                    for ( int tid = threadLB; tid < THREAD_NUM; tid++) {
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrsum_addr3( thread_r, thread_q, thread_p) == calAddrsum_addr2(r_Start, q_Start, p_Start, s_Start)) {
                                subBlkRT(cnt);
                            goto EndSample;
                        }
                    }
#ifdef DEBUG
                        if (cntStart == true) {
                            cout << "(" << calAddrsum_addr3(thread_r, thread_q, thread_p) << ", " << "(" << thread_r << ", "<< thread_q << ", "<< thread_p << ", " << "), " << cnt << ") " << endl;
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
                        int thread_r = BLIST[tid][0] + gap;
                        if (thread_r > BLIST[tid][1]) {
                            continue;
                        }
                        int thread_q = q;
                        int thread_p = p;
                        if (cntStart == true) {
                            cnt++;
#ifdef DEBUG
                            cout << "(A_addr " <<  ", " << "(" << r << ", "<< q << ", "<< p << "), " << cnt << ") " << endl;
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
int main() {
    ref_sum_addr0();
    ref_sum_addr3();
    ref_A_addr1();
    ref_A_addr0();
    ref_C4_addr0();
    ref_sum_addr1();
    ref_sum_addr2();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: doitgen */ 
