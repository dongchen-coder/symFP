
 /* Start to analysis array index
Array index info: Total number of references: 11
tmp.addr ((i * 1024) + j)
A.addr ((i * 1024) + k)
C.addr ((k * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
B.addr ((k * 1024) + j)
tmp.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
D.addr ((i * 1024) + j)
tmp.addr ((i * 1024) + k)

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %tmp
double* %A
double* %B
double* %C
double* %D
double %alpha
double %beta

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
------array access tmp.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A.addr ((i * 1024) + k)
--------array access B.addr ((k * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--------array access tmp.addr ((i * 1024) + j)
--i
--Loop Bound: (0, 1024)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 1024)
----Loop inc: (j + 1)
----Loop predicate: <
------array access D.addr ((i * 1024) + j)
------array access D.addr ((i * 1024) + j)
------k
------Loop Bound: (0, 1024)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access tmp.addr ((i * 1024) + k)
--------array access C.addr ((k * 1024) + j)
--------array access D.addr ((i * 1024) + j)
--------array access D.addr ((i * 1024) + j)

Finish analysis loops */ 
/* # of Out-most Loops: 2 */ 

 /* Start transform loop tree
	for.cond
	for.cond1
	for.cond4
	for.cond28
	for.cond31
	for.cond39
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

Finish transform loop tree */ 
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 102
------Sample number: 10485
--------Sample number: 1073741
----Sample number: 102
------Sample number: 10485
--------Sample number: 1073741
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
/* tmp_addr ((i * 1024) + j) 0 */
int calAddrtmp_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 1024) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 1024) + j) 0 */
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + j) 1 */
int calAddrtmp_addr1( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + j) 2 */
int calAddrtmp_addr2( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 0 */
int calAddrD_addr0( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 1 */
int calAddrD_addr1( int i, int j) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* tmp_addr ((i * 1024) + k) 3 */
int calAddrtmp_addr3( int i, int j, int k) {
    int result = (((i * 1024) + k)) * 8 / 64;
    return result;
}
/* C_addr ((k * 1024) + j) 0 */
int calAddrC_addr0( int i, int j, int k) {
    int result = (((k * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 2 */
int calAddrD_addr2( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
/* D_addr ((i * 1024) + j) 3 */
int calAddrD_addr3( int i, int j, int k) {
    int result = (((i * 1024) + j)) * 8 / 64;
    return result;
}
void ref_tmp_addr0() {
    cout << " ref_tmp_addr0 " << endl;
/* for (i, 0, 1024) */
/* for (j, 0, 1024) */
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
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                        cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " <<  endl;
#endif
                    }
                    while ( !thread_pool.empty()) {
                        int t_select = thread_pool[rand() % thread_pool.size()];
                        cnt++;
                        access = calAddrtmp_addr0( progress[t_select][0], progress[t_select][1]);
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
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        cnt += THREAD_NUM;
                        cnt += THREAD_NUM;
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrtmp_addr1( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrtmp_addr2( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
#ifdef DEBUG
        // cout << "Count: " << cnt << endl;
#endif
        /* Compute the chunk size. */
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    cnt += THREAD_NUM;
                    cnt += THREAD_NUM;
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrtmp_addr3( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
}
void ref_A_addr0() {
    cout << " ref_A_addr0 " << endl;
/* for (i, 0, 1024) */
/* for (j, 0, 1024) */
/* for (k, 0, 1024) */
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
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    cnt += THREAD_NUM;
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
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
                        cnt += THREAD_NUM;
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
}
void ref_C_addr0() {
    cout << " ref_C_addr0 " << endl;
/* for (i, 0, 1024) */
/* for (j, 0, 1024) */
/* for (k, 0, 1024) */
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
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    cnt += THREAD_NUM;
                    cnt += THREAD_NUM;
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        cnt += THREAD_NUM;
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrC_addr0( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
}
void ref_D_addr2() {
    cout << " ref_D_addr2 " << endl;
/* for (i, 0, 1024) */
/* for (j, 0, 1024) */
/* for (k, 0, 1024) */
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
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                        cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " <<  endl;
#endif
                    }
                    while ( !thread_pool.empty()) {
                        int t_select = thread_pool[rand() % thread_pool.size()];
                        cnt++;
                        access = calAddrD_addr0( progress[t_select][0], progress[t_select][1]);
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
                    for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                        thread_pool.push_back(*it);
#ifdef DEBUG
                        cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " <<  endl;
#endif
                    }
                    while ( !thread_pool.empty()) {
                        int t_select = thread_pool[rand() % thread_pool.size()];
                        cnt++;
                        access = calAddrD_addr1( progress[t_select][0], progress[t_select][1]);
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
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        cnt += THREAD_NUM;
                        cnt += THREAD_NUM;
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrD_addr2( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
                            thread_pool.push_back(*it);
#ifdef DEBUG
                            cout << "[" << *it << "] Iteration " << progress[*it][0] << " " << progress[*it][1] << " " << progress[*it][2] << " " <<  endl;
#endif
                        }
                        while ( !thread_pool.empty()) {
                            int t_select = thread_pool[rand() % thread_pool.size()];
                            cnt++;
                            access = calAddrD_addr3( progress[t_select][0], progress[t_select][1], progress[t_select][2]);
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
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
}
void ref_B_addr0() {
    cout << " ref_B_addr0 " << endl;
/* for (i, 0, 1024) */
/* for (j, 0, 1024) */
/* for (k, 0, 1024) */
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
        chunk_num = (1024 - 0) % (THREAD_NUM * chunk_size) == 0 ? (1024 - 0) / (THREAD_NUM * chunk_size) : (1024 - 0) / (THREAD_NUM * chunk_size) + 1;
#else
        chunk_num = 1;
        chunk_size = (1024 - 0) / THREAD_NUM;
#endif
        /* Compute the number of chunks */
        /* Generating thread local iteration space mapping code */
        for (int cid = 0; cid < chunk_num; cid++) {
            /* Computes bound express for each thread */
            for (int t = 0; t < THREAD_NUM; ++t) {
                BLIST[t][0] =  0+ (cid * THREAD_NUM + t) * chunk_size;
                BLIST[t][1] = min(0 + (cid * THREAD_NUM + t + 1) * chunk_size, 1024) - 1;
#ifdef DEBUG
                // cout << "[Thread " << t << "], " << "(" << BLIST[t][0] << ", "<< BLIST[t][1] << ")" << endl;
#endif
            }
            map<int, vector<int>> progress;
            vector<int> thread_pool;
                /* Generate the Random Interleaving process */
                vector<int> candidate_thread_pool;
                for (int tid = 0; tid < THREAD_NUM; tid++) {
                    candidate_thread_pool.push_back(tid);
                    /* init the progress vector for each thread */
                    progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0 };
                }
                while ( !candidate_thread_pool.empty()) {
                    cnt += THREAD_NUM;
                    /* Generate the Random Interleaving process */
                    vector<int> candidate_thread_pool;
                    for (int tid = 0; tid < THREAD_NUM; tid++) {
                        candidate_thread_pool.push_back(tid);
                        /* init the progress vector for each thread */
                        progress[tid] = { cid * (THREAD_NUM * chunk_size) + 0 + chunk_size * tid, 0, 0 };
                    }
                    while ( !candidate_thread_pool.empty()) {
                        cnt += THREAD_NUM;
                        for(vector<int>::iterator it = candidate_thread_pool.begin(); it != candidate_thread_pool.end(); ++it) {
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
                        cnt += THREAD_NUM;
                        cnt += THREAD_NUM;
                        for (int t_select = 0; t_select < THREAD_NUM; t_select++) {
                            /* Iteration incrementation 3 */
                            /* 2*/
                            progress[t_select][2] = ( progress[t_select][2] + 1 );
                            /* 1*/
                            progress[t_select][1] = progress[t_select][1] + (progress[t_select][2] / 1024);
                            progress[t_select][2] = progress[t_select][2] % 1024;
                            /* 0*/
                            progress[t_select][0] = progress[t_select][0] + (progress[t_select][1] / 1024);
                            progress[t_select][1] = progress[t_select][1] % 1024;
#ifdef DEBUG
                            // cout <<  "[Thread " << t_select << "] next iteration: ";
                            for (vector<int>::iterator it = progress[t_select].begin(); it != progress[t_select].end(); ++it) {
                                // cout << *it << " ";
                            }
                            // cout << endl;
#endif
                            if (progress[t_select][0] > BLIST[t_select][1]) {
                                // remove t_select from the thread pool
                                candidate_thread_pool.erase(remove(candidate_thread_pool.begin(), candidate_thread_pool.end(), t_select), candidate_thread_pool.end());
#ifdef DEBUG
                                cout << "Remove thread " << t_select << " from the candidate thread pool" << endl;
#endif
                            }
                        }
                    }
        } // end of outer for - cid loops
}
int main() {
    ref_tmp_addr0();
    LAT.clear();
    ref_A_addr0();
    LAT.clear();
    ref_C_addr0();
    LAT.clear();
    ref_D_addr2();
    LAT.clear();
    ref_B_addr0();
    LAT.clear();
    rtDump();
    RTtoMR_AET();
    dumpMR();
    return 0;
}
 /* Analyze function: mm2 */ 
