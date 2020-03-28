
 /* Start to analysis array index
Array index info: Total number of references: 1
A.addr i

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B
double %alpha
double %beta

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 32)
--Loop inc: (i + 1)
--Loop predicate: <
----r
----Loop Bound: (0, 16)
----Loop inc: (r + 1)
----Loop predicate: <
------j
------Loop Bound: (0, 32)
------Loop inc: (j + 8)
------Loop predicate: <
--------array access A.addr i

Finish analysis loops */ 
/* # of Out-most Loops: 1 */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 32
------Sample number: 512
--------Sample number: 2048
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
/* A_addr0	2048 */
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <cmath>
#ifdef PAPI_TIMER
#  include <chrono>
#endif
using namespace std;
#ifdef PAPI_TIMER
using namespace  chrono;
#endif
 map<uint64_t, double> RT;
 map<uint64_t, double> MR;
void rtHistoCal( map<uint64_t, double> &rth, uint64_t rt, double val ) {
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
    }
    return;
}
void subBlkRT(map<uint64_t, double> &rth, int rt, double cnt) {
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
                rtHistoCal(rth, b - diff, cnt);
                break;
            }
        }
    }
    else {
        rtHistoCal(rth, pow(2, msb-1), cnt);
    }
    return;
}
void RTtoMR_AET() {
     map<uint64_t, double> P;
    double total_num_RT = 0;
    uint64_t max_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
        total_num_RT += it->second;
        if (max_RT < it->first) {
            max_RT = it->first;
        }
    }
    double accumulate_num_RT = 0;
    for ( map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {
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
        cout << it->first << ", " << it->second << "\n";
    }
    return;
}
void dumpMR() {
    cout << "miss ratio" << endl;
     map<uint64_t, double>::iterator it1 = MR.begin();
     map<uint64_t, double>::iterator it2 = MR.begin();
    while(it1 != MR.end()) {
        while(1) {
             map<uint64_t, double>::iterator it3 = it2;
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
        cout << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
/* Array A_addr	i */ 
/* A_addr i 0 */
int calAddrA_addr0( int i, int r, int j) {
    int result = (i) * 8 / 64;
    return result;
}
int getChunkID(uint64_t i) {
    return floor(i / (CHUNK_SIZE * THREAD_NUM));
}
int getThreadID(uint64_t i) {
    return i / CHUNK_SIZE - floor(i / (CHUNK_SIZE * THREAD_NUM))*THREAD_NUM ;
}
int getThreadLocalPos(uint64_t i) {
    return i % CHUNK_SIZE;
}
int search_src_candidate_neighbor(uint64_t i, uint64_t r, uint64_t j) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    int c_end = c_start + (THREAD_NUM - 1) * CHUNK_SIZE;
    for (int c = i + CHUNK_SIZE; c <= c_end; c = c + CHUNK_SIZE) {
        if (calAddrA_addr0(i, r, j) == calAddrA_addr0(c, r, j)) {
            return getThreadID(c); 
        }
    }
    return -1;
}
int search_sink_candidate_neighbor(uint64_t i, uint64_t r, uint64_t j) {
    int c_start = i % CHUNK_SIZE + i / (CHUNK_SIZE * THREAD_NUM) * THREAD_NUM * CHUNK_SIZE;
    for (int c = c_start; c <= i; c = c + CHUNK_SIZE) {
        if (calAddrA_addr0(i, r, j) == calAddrA_addr0(c, r, j)) {
            return getThreadID(c); 
        }
    }
    return -1;
}
uint64_t parallel_predict(uint64_t i_src, uint64_t r_src, uint64_t j_src, uint64_t i_sink, uint64_t r_sink, uint64_t j_sink, uint64_t rt, uint64_t lsrc, uint64_t lsink, bool is_normal_ref) {
    uint64_t parallel_rt = rt;
    int tsrc = getThreadID(i_src);
    int tsink = getThreadID(i_sink);
    int dT = tsink - tsrc;
    if (!is_normal_ref) {
#ifdef DEBUG
        cout << "Neighboring Effect" << endl;
#endif
        int tneighbor_src = search_src_candidate_neighbor(i_src, r_src, j_src);
        int tneighbor_sink = search_sink_candidate_neighbor(i_sink, r_sink, j_sink);
        if (tneighbor_src >= 0) {
#ifdef DEBUG
        cout << "Find sink in src neighbor at " << tneighbor_src << endl;
#endif
            return tneighbor_src - tsrc;
        } else if (tneighbor_sink >= 0) {
#ifdef DEBUG
        cout << "Find sink in sink neighbor at " << tneighbor_sink << endl;
#endif
            return rt * THREAD_NUM + tneighbor_sink - tsink;
        }
    }
    /* intra chunk reuse */
    if (getChunkID(i_src) == getChunkID(i_sink)) {
        /* same thread -- scaling effect */
        if (dT == 0) {
#ifdef DEBUG
            cout << "Scaling Effect" << endl;
#endif
            parallel_rt = rt * THREAD_NUM;
        } else if (getThreadLocalPos(i_src) <= getThreadLocalPos(i_sink)) { // src-sink order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("NORMAL ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Src-Sink Order Folding Effect" << endl;
#endif
            parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + abs(dT);
        } else { // sink-src order
            if ((rt * THREAD_NUM - CHUNK_SIZE * lsrc * THREAD_NUM * dT + dT) < 0) { printf("REVERSE ORDER NEGATIVE PRI\n"); }
#ifdef DEBUG
            cout << "Sink-Src Order Folding Effect" << endl;
#endif
            parallel_rt = CHUNK_SIZE * lsrc * THREAD_NUM * dT - (rt * THREAD_NUM) - abs(dT);
        }
    } else { // inter chunk reuse 
#ifdef DEBUG
            cout << "Inter Chunk Reuse" << endl;
#endif
            parallel_rt = rt * THREAD_NUM - CHUNK_SIZE * THREAD_NUM * (lsrc*(THREAD_NUM - tsrc) + lsink * tsink) + CHUNK_SIZE * THREAD_NUM * lsink + dT;
    }
    return parallel_rt;
}

void ref_A_addr0() {
    /* Generating sampling loop */
    set<string> record;
    for ( int s = 0; s < 2048;) {
SAMPLE:
        int i_Start = rand() % (32 - 0) + 0;
        if (i_Start % 1 != 0) goto SAMPLE; 
        if ( (16 - 0) == 0) goto SAMPLE;
        int r_Start = rand() % (16 - 0) + 0;
        if (r_Start % 1 != 0) goto SAMPLE; 
        if ( (32 - 0) == 0) goto SAMPLE;
        int j_Start = rand() % (32 - 0) + 0;
        if (j_Start % 8 != 0) goto SAMPLE; 
        string idx_string =  to_string(i_Start) + "_" +  to_string(r_Start) + "_" +  to_string(j_Start) + "_" ;
        if ( record.find(idx_string) != record.end() ) goto SAMPLE;
        record.insert( idx_string );
        uint64_t cnt = 0;
        bool cntStart = false;

        /* Generating reuse search code */

        {
        int iLB0 = i_Start;
        for ( int i = iLB0; i < 32; i=(i + 1)) {
            {
            int rLB1 = 0;
            if ( i == i_Start ) {
                rLB1 = r_Start;
            }
            for ( int r = rLB1; r < 16; r=(r + 1)) {
                {
                int jLB2 = 0;
                if ( i == i_Start && r == r_Start ) {
                    jLB2 = j_Start;
                }
                for ( int j = jLB2; j < 32; j=(j + 8)) {
                    if (cntStart == true) {
                        cnt++;
                        if ( calAddrA_addr0( i, r, j) == calAddrA_addr0(i_Start, r_Start, j_Start)) {
                                uint64_t parallel_rt = parallel_predict((i_Start -0), r_Start, j_Start, i, r, j, cnt, 2048, 2048, false);
                                rtHistoCal(RT, parallel_rt, 1.0);
#ifdef DEBUG
                                    cout << "[" << parallel_rt << "] (" << i_Start<< ", " << r_Start<< ", " << j_Start<< ") -- (" << i<< ", " << r<< ", " << j<< ") " << endl;
#endif
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
#ifdef PAPI_TIMER
    // Get starting timepoint
    auto start = high_resolution_clock::now();
#endif
    ref_A_addr0();
    rtDump();
    RTtoMR_AET();
    dumpMR();
#ifdef PAPI_TIMER
    // Get ending timepoint
    auto stop = high_resolution_clock::now(); 
    // Get duration. Substart timepoints to
    // get durarion. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);
     cout << "Time taken by SPS:  " << duration.count() << endl; 
#endif
    return 0;
}
 /* Analyze function: test_kernel */ 
