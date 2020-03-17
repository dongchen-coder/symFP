
 /* Start to analysis array index
Array index info: Total number of references: 8
C.addr ((i * 128) + j)
C.addr ((i * 128) + j)
B.addr ((j * 128) + k)
A.addr ((i * 128) + k)
A.addr ((j * 128) + k)
B.addr ((i * 128) + k)
C.addr ((i * 128) + j)
C.addr ((i * 128) + j)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %A
double* %B
double* %C
double %alpha
double %beta

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------array access C.addr ((i * 128) + j)
------array access C.addr ((i * 128) + j)
--i
--Loop Bound: (0, 128)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 128)
----Loop inc: (j + 1)
----Loop predicate: <
------k
------Loop Bound: (0, 128)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((j * 128) + k)
--------array access B.addr ((i * 128) + k)
--------array access B.addr ((j * 128) + k)
--------array access A.addr ((i * 128) + k)
--------array access C.addr ((i * 128) + j)
--------array access C.addr ((i * 128) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 2 */ 

 /* Start transform loop tree
	for.cond
	for.cond1
	for.cond8
	for.cond11
	for.cond14
----------------
--|  LoopNode  |
----------------
------------------
----|  LoopNode  |
------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
--------------------
------| ThreadNode |
--------------------
----------------------
--------| AccessNode |
----------------------
----------------
--|  LoopNode  |
----------------
------------------
----|  LoopNode  |
------------------
--------------------
------|  LoopNode  |
--------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------
----------------------
--------| ThreadNode |
----------------------
------------------------
----------| AccessNode |
------------------------

Finish transform loop tree */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 12
------Sample number: 163
----Sample number: 12
------Sample number: 163
--------Sample number: 2097
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cmath>
#ifndef THREAD_NUM
#    define THREAD_NUM   4
#endif
#ifndef BIN_SIZE
#    define BIN_SIZE   4
#endif
#ifndef CHUNK_SIZE
#    define CHUNK_SIZE   4
#endif
using namespace std;
std::map<uint64_t, tuple<uint64_t, int>> LAT;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( map<uint64_t, double> &rth, int rt, int val ) {
    if ( val <= 0) {
;        return;
    }
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
    }
    return;
}
void subBlkRT(map<uint64_t, double> &rth, int rt) {
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
                rtHistoCal(rth, b - diff, 1);
                break;
            }
        }
    }
    else {
        rtHistoCal(rth, pow(2, msb-1), 1);
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
        cout << it->first << ", " << it->second << "\n";
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
        cout << it1->first << ", " << it1->second << endl;
        if (it1 != it2) {
            cout << it2->first << ", " << it2->second << endl;
        }
        it1 = ++it2;
        it2 = it1;
    }
    return;
}
/* C_addr ((i * 128) + j) 0 */
int calAddrC_addr0( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 128) + j) 1 */
int calAddrC_addr1( int i, int j) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* A_addr ((j * 128) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((j * 128) + k)) * 8 / 64;
    return result;
}
/* B_addr ((i * 128) + k) 0 */
int calAddrB_addr0( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* B_addr ((j * 128) + k) 1 */
int calAddrB_addr1( int i, int j, int k) {
    int result = (((j * 128) + k)) * 8 / 64;
    return result;
}
/* A_addr ((i * 128) + k) 1 */
int calAddrA_addr1( int i, int j, int k) {
    int result = (((i * 128) + k)) * 8 / 64;
    return result;
}
/* C_addr ((i * 128) + j) 2 */
int calAddrC_addr2( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 128) + j) 3 */
int calAddrC_addr3( int i, int j, int k) {
    int result = (((i * 128) + j)) * 8 / 64;
    return result;
}
void ref_C_addr0() {
    cout << " ref_C_addr0 " << endl;
/* for (i, 0, 128) */
/* for (j, 0, 128) */
    uint64_t cnt = 0;

    /* Variable used to compute thread-local iteration space (out-most-loops) */
    auto BLIST = new int[THREAD_NUM][2];
    int t_Start = 0;
    /* Generating reuse search code */
    /* Sampled IDVs 2  */
    /* Sampled IDV: i  */
    /* Sampled IDV: j  */
    /* Vector that contains the interleaved iteration, avoid duplicate declaration */
    vector<vector<int>> nv(THREAD_NUM);
    int chunk_size, chunk_num;
    uint64_t access;
#ifdef DEBUG
    // cout << "Count: " << cnt << endl;
#endif
    /* Compute the chunk size. */
#ifdef CHUNK_SIZE
    chunk_size = CHUNK_SIZE;
    chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
    chunk_num = 1;
    chunk_size = (128 - 0) / THREAD_NUM;
#endif
    /* Compute the number of chunks */
    /* Generating thread local iteration space mapping code */
    for (int cid = 0; cid < chunk_num; cid++) {
        /* Computes bound express for each thread */
        for (int t = 0; t < THREAD_NUM; ++t) {
            BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
            BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
            // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
        }
        vector<int> thread_pool;
        map<int, vector<int>> progress;
        /* Generate the Random Interleaving process */
        vector<int> candidate_thread_pool_2;
        for (int tid = 0; tid < THREAD_NUM; tid++) {
            candidate_thread_pool_2.push_back(tid);
            /* init the progress vector for each thread (1) */
            progress[tid] = {cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
        }
        while ( !candidate_thread_pool_2.empty()) {
            for(vector<int>::iterator it = candidate_thread_pool_2.begin(); it != candidate_thread_pool_2.end(); ++it) {
                    thread_pool.push_back(*it);
#ifdef DEBUG
                cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " <<  endl;
#endif
            }
            while ( !thread_pool.empty()) {
                int t_select = thread_pool[rand() % thread_pool.size()];
                cnt++;
                access = calAddrC_addr0( progress[t_select][0], progress[t_select][1]);
                if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                    cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< ")" << endl;
#endif
                    subBlkRT(RT, cnt - get<0>(LAT[access]));
                }
                LAT[access] = make_tuple(cnt, cid);
                thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
            }
            for(vector<int>::iterator it = candidate_thread_pool_2.begin(); it != candidate_thread_pool_2.end(); ++it) {
                    thread_pool.push_back(*it);
#ifdef DEBUG
                cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " <<  endl;
#endif
            }
            while ( !thread_pool.empty()) {
                int t_select = thread_pool[rand() % thread_pool.size()];
                cnt++;
                access = calAddrC_addr1( progress[t_select][0], progress[t_select][1]);
                if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                    cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< ")" << endl;
#endif
                    subBlkRT(RT, cnt - get<0>(LAT[access]));
                }
                LAT[access] = make_tuple(cnt, cid);
                thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
            }
            /* hasTNLoop: True */
            /* containsFirstTN: False */
            /* !containsFirstTN && hasTNLoop: True */
            for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                /* Iteration incrementation 1 */
                progress[t_select][1] = progress[t_select][1] + 1; 
                progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 128);
                progress[t_select][1] = progress[t_select][1] % 128;
#ifdef DEBUG
                // cout <<  "[Thread " << t_select << "] next iteration: ";
                for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                    // cout << *it << " ";
                }
                // cout << endl;
#endif
                if (progress[t_select][0] > BLIST[t_select][1]) {
                    // remove t_select from the thread pool
                    candidate_thread_pool_2.erase(remove(candidate_thread_pool_2.begin(), candidate_thread_pool_2.end(), t_select), candidate_thread_pool_2.end());
#ifdef DEBUG
                    cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                }
            }
        } // end of while loop
    } // end of outer for - cid loops
#ifdef DEBUG
        // cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (128 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            vector<int> thread_pool;
            map<int, vector<int>> progress;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool_3;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool_3.push_back(tid);
                    /* init the progress vector for each thread (2) */
                    progress[tid] = {cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                }
                while ( !candidate_thread_pool_3.empty()) {
                    cnt += THREAD_NUM;
                    cnt += THREAD_NUM;
                    cnt += THREAD_NUM;
                    cnt += THREAD_NUM;
                    for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                        cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                    }
                    while ( !thread_pool.empty()) {
                        int t_select = thread_pool[rand() % thread_pool.size()];
                        cnt++;
                        access = calAddrC_addr2( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                        if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                            cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                            subBlkRT(RT, cnt - get<0>(LAT[access]));
                        }
                        LAT[access] = make_tuple(cnt, cid);
                        thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                        // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                    }
                    for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                        cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                    }
                    while ( !thread_pool.empty()) {
                        int t_select = thread_pool[rand() % thread_pool.size()];
                        cnt++;
                        access = calAddrC_addr3( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                        if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                            cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                            subBlkRT(RT, cnt - get<0>(LAT[access]));
                        }
                        LAT[access] = make_tuple(cnt, cid);
                        thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                        // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                    }
                    /* hasTNLoop: True */
                    /* containsFirstTN: False */
                    /* !containsFirstTN && hasTNLoop: True */
                    for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                        /* Iteration incrementation 2 */
                        progress[t_select][2] = progress[t_select][2] + 1; 
                        progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 128);
                        progress[t_select][2] = progress[t_select][2] % 128;
                        progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 128);
                        progress[t_select][1] = progress[t_select][1] % 128;
#ifdef DEBUG
                        // cout <<  "[Thread " << t_select << "] next iteration: ";
                        for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                            // cout << *it << " ";
                        }
                        // cout << endl;
#endif
                        if (progress[t_select][0] > BLIST[t_select][1]) {
                            // remove t_select from the thread pool
                            candidate_thread_pool_3.erase(remove(candidate_thread_pool_3.begin(), candidate_thread_pool_3.end(), t_select), candidate_thread_pool_3.end());
#ifdef DEBUG
                            cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                        }
                    }
                } // end of while loop
        } // end of outer for - cid loops
}
void ref_B_addr1() {
    cout << " ref_B_addr1 " << endl;
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    uint64_t cnt = 0;

    /* Variable used to compute thread-local iteration space (out-most-loops) */
    auto BLIST = new int[THREAD_NUM][2];
    int t_Start = 0;
    /* Generating reuse search code */
    /* Sampled IDVs 3  */
    /* Sampled IDV: i  */
    /* Sampled IDV: j  */
    /* Sampled IDV: k  */
    /* Vector that contains the interleaved iteration, avoid duplicate declaration */
    vector<vector<int>> nv(THREAD_NUM);
    int chunk_size, chunk_num;
    uint64_t access;
#ifdef DEBUG
    // cout << "Count: " << cnt << endl;
#endif
    /* Compute the chunk size. */
#ifdef CHUNK_SIZE
    chunk_size = CHUNK_SIZE;
    chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
    chunk_num = 1;
    chunk_size = (128 - 0) / THREAD_NUM;
#endif
    /* Compute the number of chunks */
    /* Generating thread local iteration space mapping code */
    for (int cid = 0; cid < chunk_num; cid++) {
        /* Computes bound express for each thread */
        for (int t = 0; t < THREAD_NUM; ++t) {
            BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
            BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
            // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
        }
        vector<int> thread_pool;
        map<int, vector<int>> progress;
            /* Generate the Random Interleaving process */
            vector<int> candidate_thread_pool_3;
            for (int tid = 0; tid < THREAD_NUM; tid++) {
                candidate_thread_pool_3.push_back(tid);
                /* init the progress vector for each thread (2) */
                progress[tid] = {cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
            }
            while ( !candidate_thread_pool_3.empty()) {
                cnt += THREAD_NUM;
                for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                    cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                }
                while ( !thread_pool.empty()) {
                    int t_select = thread_pool[rand() % thread_pool.size()];
                    cnt++;
                    access = calAddrB_addr0( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                    if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                        cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                        subBlkRT(RT, cnt - get<0>(LAT[access]));
                    }
                    LAT[access] = make_tuple(cnt, cid);
                    thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                    // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                }
                for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                    cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                }
                while ( !thread_pool.empty()) {
                    int t_select = thread_pool[rand() % thread_pool.size()];
                    cnt++;
                    access = calAddrB_addr1( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                    if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                        cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                        subBlkRT(RT, cnt - get<0>(LAT[access]));
                    }
                    LAT[access] = make_tuple(cnt, cid);
                    thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                    // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                }
                cnt += THREAD_NUM;
                cnt += THREAD_NUM;
                cnt += THREAD_NUM;
                /* hasTNLoop: True */
                /* containsFirstTN: False */
                /* !containsFirstTN && hasTNLoop: True */
                for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                    /* Iteration incrementation 2 */
                    progress[t_select][2] = progress[t_select][2] + 1; 
                    progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 128);
                    progress[t_select][2] = progress[t_select][2] % 128;
                    progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 128);
                    progress[t_select][1] = progress[t_select][1] % 128;
#ifdef DEBUG
                    // cout <<  "[Thread " << t_select << "] next iteration: ";
                    for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                        // cout << *it << " ";
                    }
                    // cout << endl;
#endif
                    if (progress[t_select][0] > BLIST[t_select][1]) {
                        // remove t_select from the thread pool
                        candidate_thread_pool_3.erase(remove(candidate_thread_pool_3.begin(), candidate_thread_pool_3.end(), t_select), candidate_thread_pool_3.end());
#ifdef DEBUG
                        cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                    }
                }
            } // end of while loop
    } // end of outer for - cid loops
}
void ref_A_addr1() {
    cout << " ref_A_addr1 " << endl;
/* for (i, 0, 128) */
/* for (j, 0, 128) */
/* for (k, 0, 128) */
    uint64_t cnt = 0;

    /* Variable used to compute thread-local iteration space (out-most-loops) */
    auto BLIST = new int[THREAD_NUM][2];
    int t_Start = 0;
    /* Generating reuse search code */
    /* Sampled IDVs 3  */
    /* Sampled IDV: i  */
    /* Sampled IDV: j  */
    /* Sampled IDV: k  */
    /* Vector that contains the interleaved iteration, avoid duplicate declaration */
    vector<vector<int>> nv(THREAD_NUM);
    int chunk_size, chunk_num;
    uint64_t access;
#ifdef DEBUG
    // cout << "Count: " << cnt << endl;
#endif
    /* Compute the chunk size. */
#ifdef CHUNK_SIZE
    chunk_size = CHUNK_SIZE;
    chunk_num = (128 - 0) % (THREAD_NUM * chunk_size) == 0 ? (128 - 0) / (THREAD_NUM * chunk_size) : (128 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
    chunk_num = 1;
    chunk_size = (128 - 0) / THREAD_NUM;
#endif
    /* Compute the number of chunks */
    /* Generating thread local iteration space mapping code */
    for (int cid = 0; cid < chunk_num; cid++) {
        /* Computes bound express for each thread */
        for (int t = 0; t < THREAD_NUM; ++t) {
            BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
            BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 128) - 1;
#ifdef DEBUG
            // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
        }
        vector<int> thread_pool;
        map<int, vector<int>> progress;
            /* Generate the Random Interleaving process */
            vector<int> candidate_thread_pool_3;
            for (int tid = 0; tid < THREAD_NUM; tid++) {
                candidate_thread_pool_3.push_back(tid);
                /* init the progress vector for each thread (2) */
                progress[tid] = {cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
            }
            while ( !candidate_thread_pool_3.empty()) {
                for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                    cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                }
                while ( !thread_pool.empty()) {
                    int t_select = thread_pool[rand() % thread_pool.size()];
                    cnt++;
                    access = calAddrA_addr0( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                    if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                        cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                        subBlkRT(RT, cnt - get<0>(LAT[access]));
                    }
                    LAT[access] = make_tuple(cnt, cid);
                    thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                    // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                }
                cnt += THREAD_NUM;
                cnt += THREAD_NUM;
                for(vector<int>::iterator it = candidate_thread_pool_3.begin(); it != candidate_thread_pool_3.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                    cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                }
                while ( !thread_pool.empty()) {
                    int t_select = thread_pool[rand() % thread_pool.size()];
                    cnt++;
                    access = calAddrA_addr1( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
                    if (LAT.find(access) != LAT.end()) {
#ifdef DEBUG
                        cout << "[REUSE of Addr " << access << "] " << cnt - get<0>(LAT[access]) << " find @(" << progress[t_select][0] << " "<< progress[t_select][1] << " "<< progress[t_select][2] << " "<< ")" << endl;
#endif
                        subBlkRT(RT, cnt - get<0>(LAT[access]));
                    }
                    LAT[access] = make_tuple(cnt, cid);
                    thread_pool.erase(remove(thread_pool.begin(), thread_pool.end(), t_select), thread_pool.end());
#ifdef DEBUG
                    // cout << "Remove thread " << t_select << " from the pool" << endl;
#endif
                }
                cnt += THREAD_NUM;
                cnt += THREAD_NUM;
                /* hasTNLoop: True */
                /* containsFirstTN: False */
                /* !containsFirstTN && hasTNLoop: True */
                for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                    /* Iteration incrementation 2 */
                    progress[t_select][2] = progress[t_select][2] + 1; 
                    progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 128);
                    progress[t_select][2] = progress[t_select][2] % 128;
                    progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 128);
                    progress[t_select][1] = progress[t_select][1] % 128;
#ifdef DEBUG
                    // cout <<  "[Thread " << t_select << "] next iteration: ";
                    for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                        // cout << *it << " ";
                    }
                    // cout << endl;
#endif
                    if (progress[t_select][0] > BLIST[t_select][1]) {
                        // remove t_select from the thread pool
                        candidate_thread_pool_3.erase(remove(candidate_thread_pool_3.begin(), candidate_thread_pool_3.end(), t_select), candidate_thread_pool_3.end());
#ifdef DEBUG
                        cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                    }
                }
            } // end of while loop
    } // end of outer for - cid loops
}
int main() {
    ref_C_addr0();
    LAT.clear();
    ref_B_addr1();
    LAT.clear();
    ref_A_addr1();
    LAT.clear();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: syr2k */ 
