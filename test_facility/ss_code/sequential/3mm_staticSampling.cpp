
 /* Start to analysis array index
Array index info: Total number of references: 15
E.addr ((i * 2048) + j)
E.addr ((i * 2048) + j)
E.addr ((i * 2048) + j)
C.addr ((i * 2048) + k)
D.addr ((k * 2048) + j)
A.addr ((i * 2048) + k)
B.addr ((k * 2048) + j)
G.addr ((i * 2048) + j)
F.addr ((i * 2048) + j)
F.addr ((i * 2048) + j)
F.addr ((i * 2048) + j)
E.addr ((i * 2048) + k)
F.addr ((k * 2048) + j)
G.addr ((i * 2048) + j)
G.addr ((i * 2048) + j)
BC Array cost info: Total number of arrays: 7
A.addr 2
B.addr 2
C.addr 2
D.addr 2
E.addr 8
F.addr 8
G.addr 6
BC Average cost: 4.285714e+00

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %ni
i32 %nj
i32 %nk
i32 %nl
i32 %nm
double* %E
double* %A
double* %B
double* %F
double* %C
double* %D
double* %G

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
Array F_addr 
Array B_addr 
Array D_addr 
Array G_addr 
Array A_addr 
Array C_addr 
Array E_addr 
--i
--Loop Bound: (0, 2048)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 2048)
----Loop inc: (j + 1)
----Loop predicate: <
------array access E_addr ((i * 2048) + j)
------k
------Loop Bound: (0, 2048)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access A_addr ((i * 2048) + k)
--------array access B_addr ((k * 2048) + j)
--------array access E_addr ((i * 2048) + j)
--------array access E_addr ((i * 2048) + j)
--i
--Loop Bound: (0, 2048)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 2048)
----Loop inc: (j + 1)
----Loop predicate: <
------array access F_addr ((i * 2048) + j)
------k
------Loop Bound: (0, 2048)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access C_addr ((i * 2048) + k)
--------array access D_addr ((k * 2048) + j)
--------array access F_addr ((i * 2048) + j)
--------array access F_addr ((i * 2048) + j)
--i
--Loop Bound: (0, 2048)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 2048)
----Loop inc: (j + 1)
----Loop predicate: <
------array access G_addr ((i * 2048) + j)
------k
------Loop Bound: (0, 2048)
------Loop inc: (k + 1)
------Loop predicate: <
--------array access E_addr ((i * 2048) + k)
--------array access F_addr ((k * 2048) + j)
--------array access G_addr ((i * 2048) + j)
--------array access G_addr ((i * 2048) + j)

Finish analysis loops */ 
/* Start IV Dependence Analysis Finish to analyze IV Dependence *//* Start to analysis the access graphG_addr -> E_addr
E_addr -> F_addr
F_addr -> G_addr
G_addr -> G_addr
G_addr -> E_addr
G_addr -> G_addr
G_addr -> G_addr

F_addr -> C_addr
C_addr -> D_addr
D_addr -> F_addr
F_addr -> F_addr
F_addr -> C_addr
F_addr -> F_addr
F_addr -> F_addr

E_addr -> A_addr
A_addr -> B_addr
B_addr -> E_addr
E_addr -> E_addr
E_addr -> A_addr
E_addr -> E_addr
E_addr -> E_addr

*/
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 204
------Sample number: 41943
--------Sample number: 8589934
----Sample number: 204
------Sample number: 41943
--------Sample number: 8589934
----Sample number: 204
------Sample number: 41943
--------Sample number: 8589934
 End of sample analysis */
 // Start to generating Static Sampling Code (reference based)
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <time.h>
#include <cassert>
#ifdef PAPI_TIMER
#  include "papi_timer.h"
#endif
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
typedef pair<int, int> Chunk;
class ChunkEngine {
    int lb = 0;
    int ub = 0;
    int chunk_size = 0;
    int trip = 0;
    int avail_chunk = 0;
public:
    ChunkEngine() {} 
    ChunkEngine(int chunk_size, int trip) {
        assert(chunk_size <= trip);
        this->chunk_size = chunk_size;
        this->trip = trip;
        this->avail_chunk = (trip / chunk_size + (trip % chunk_size));
        this->lb = 0;
        this->ub = chunk_size - 1;
    }
    string getCurrentChunkRange() {
        return "[" + to_string(this->lb) + ", " + to_string(this->ub) + "]";
    }
    bool hasNextChunk() {
        return this->avail_chunk > 0;
    }
    Chunk getNextChunk(int tid) {
        // assign the current lb, ub to thread tid and update the next chunk
        Chunk curr = make_pair(this->lb, this->ub);
        this->lb = this->ub + 1;
        this->ub = (this->lb + chunk_size - 1) <= this->trip ? (this->lb + chunk_size - 1) : this->trip;
        this->avail_chunk -= 1;
        return curr;
    }
};
class Progress {
public:
    string ref;
    Chunk chunk;
    vector<int> iteration;
    Progress() { }
    Progress(string ref, vector<int> iteration, Chunk c) {
        this->ref = ref;
        this->iteration = iteration;
        this->chunk = c;
    }
    string getIteration() {
        string ret = "(";
        for (int i = 0; i < this->iteration.size(); i++) {
            ret += to_string(this->iteration[i]);
            if (i != this->iteration.size() - 1)
                ret += ",";
            }
        ret += ")";
        return ret;
    }
    string getReference() {
        return this->ref;
    }
    void increment(string ref, vector<int> iteration) {
        this->ref = ref;
        this->iteration = iteration;
    }
    void increment(string ref) {
        this->ref = ref;
    }
    bool isInBound() {
        assert(this->iteration[0] >= chunk.first);
        return this->iteration[0] <= chunk.second;
    }
};
std::unordered_map<uint64_t, uint64_t> LAT_E_addr;
std::unordered_map<uint64_t, uint64_t> LAT_C_addr;
std::unordered_map<uint64_t, uint64_t> LAT_D_addr;
std::unordered_map<uint64_t, uint64_t> LAT_A_addr;
std::unordered_map<uint64_t, uint64_t> LAT_B_addr;
std::unordered_map<uint64_t, uint64_t> LAT_G_addr;
std::unordered_map<uint64_t, uint64_t> LAT_F_addr;
std::map<uint64_t, double> RT;
std::map<uint64_t, double> MR;
void rtHistoCal( map<uint64_t, double> &rth, int rt, int val ) {
    if ( val <= 0) {
        return;
    }
    if (rth.find(rt) == rth.end()) { 
        rth[rt] = val;
    } else {
        rth[rt] += val;
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
    cout << fixed << setprecision(3);
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
/* E_addr ((i * 2048) + j) 0 */
int calAddrE_addr0( int i, int j) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* A_addr ((i * 2048) + k) 0 */
int calAddrA_addr0( int i, int j, int k) {
    int result = (((i * 2048) + k)) * 8 / 64;
    return result;
}
/* B_addr ((k * 2048) + j) 0 */
int calAddrB_addr0( int i, int j, int k) {
    int result = (((k * 2048) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 2048) + j) 1 */
int calAddrE_addr1( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 2048) + j) 2 */
int calAddrE_addr2( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 2048) + j) 0 */
int calAddrF_addr0( int i, int j) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* C_addr ((i * 2048) + k) 0 */
int calAddrC_addr0( int i, int j, int k) {
    int result = (((i * 2048) + k)) * 8 / 64;
    return result;
}
/* D_addr ((k * 2048) + j) 0 */
int calAddrD_addr0( int i, int j, int k) {
    int result = (((k * 2048) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 2048) + j) 1 */
int calAddrF_addr1( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* F_addr ((i * 2048) + j) 2 */
int calAddrF_addr2( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 2048) + j) 0 */
int calAddrG_addr0( int i, int j) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* E_addr ((i * 2048) + k) 3 */
int calAddrE_addr3( int i, int j, int k) {
    int result = (((i * 2048) + k)) * 8 / 64;
    return result;
}
/* F_addr ((k * 2048) + j) 3 */
int calAddrF_addr3( int i, int j, int k) {
    int result = (((k * 2048) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 2048) + j) 1 */
int calAddrG_addr1( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
/* G_addr ((i * 2048) + j) 2 */
int calAddrG_addr2( int i, int j, int k) {
    int result = (((i * 2048) + j)) * 8 / 64;
    return result;
}
void interleaving() {
    uint64_t cnt = 0;
    int tid_to_run, chunk_size;
    vector<int> candidate_thread_pool, threads_to_exec;
    auto randgen = []() { return rand() % THREAD_NUM; };
    ChunkEngine engine;
    unordered_map<int, Progress> progress;
#ifdef CHUNK_SIZE
        chunk_size = CHUNK_SIZE;
#else
        chunk_size = 1;
#endif
    /* USEING STATIC SCHEUDLING */
    engine = ChunkEngine(chunk_size, (2048 - 0));
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        candidate_thread_pool.emplace_back(tid_to_run);
    }
    while(true) {
        if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {
            while(!candidate_thread_pool.empty()) {
                tid_to_run = *(candidate_thread_pool.begin());
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] Assign chunk " << engine.getCurrentChunkRange() << endl;
#endif
                Chunk c = engine.getNextChunk(tid_to_run);
                Progress p("E_addr0", {c.first, 0}, c);
                progress[tid_to_run] = p;
                candidate_thread_pool.erase(candidate_thread_pool.begin());
                threads_to_exec.emplace_back(tid_to_run);
            } /* end of progress assignment */
        } /* end of chunk availability check */
        /* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
            if (!progress[tid_to_run].isInBound()) {
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] " << progress[tid_to_run].iteration[0] << " > " << progress[tid_to_run].chunk.second << endl;
#endif
                continue;
            }
#ifdef SIMULATOR_DEBUG
            cout << "[" << tid_to_run << "] Iterate " << progress[tid_to_run].ref << " at " << progress[tid_to_run].getIteration() << endl;
#endif
            if (progress[tid_to_run].ref == "E_addr0") {
                cnt++;
                int access = calAddrE_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_E_addr.find(access) != LAT_E_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_E_addr[access], 1.0);
                    }
                }
                LAT_E_addr[access] = cnt;
/* E_addr0 and A_addr0 are not in the same loop level */
/*   %k = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("A_addr0");
                continue;
            } /* end of check to E_addr0 */
            if (progress[tid_to_run].ref == "A_addr0") {
                cnt++;
                int access = calAddrA_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_A_addr.find(access) != LAT_A_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_A_addr[access], 1.0);
                    }
                }
                LAT_A_addr[access] = cnt;
                progress[tid_to_run].increment("B_addr0");
                continue;
            } /* end of check to A_addr0 */
            if (progress[tid_to_run].ref == "B_addr0") {
                cnt++;
                int access = calAddrB_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_B_addr.find(access) != LAT_B_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_B_addr[access], 1.0);
                    }
                }
                LAT_B_addr[access] = cnt;
                progress[tid_to_run].increment("E_addr1");
                continue;
            } /* end of check to B_addr0 */
            if (progress[tid_to_run].ref == "E_addr1") {
                cnt++;
                int access = calAddrE_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_E_addr.find(access) != LAT_E_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_E_addr[access], 1.0);
                    }
                }
                LAT_E_addr[access] = cnt;
                progress[tid_to_run].increment("E_addr2");
                continue;
            } /* end of check to E_addr1 */
            if (progress[tid_to_run].ref == "E_addr2") {
                cnt++;
                int access = calAddrE_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_E_addr.find(access) != LAT_E_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_E_addr[access], 1.0);
                    }
                }
                LAT_E_addr[access] = cnt;
/* E_addr2 is the last access node in a loop */
/*   %k = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[2] + 1 < 2048) {
                    progress[tid_to_run].iteration[2] += 1;
                    progress[tid_to_run].increment("A_addr0");
                    continue;
                } /* end of check to E_addr2 */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 2048) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("E_addr0");
                    continue;
                } /* end of check to E_addr2 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].iteration.pop_back();
/*   %j = alloca i32, align 4 has no depend parent */
                    progress[tid_to_run].iteration.emplace_back(0);
                    progress[tid_to_run].increment("E_addr0");
                    continue;
                } /* end of check to E_addr2 */
                if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {
                    candidate_thread_pool.emplace_back(tid_to_run);
                }
            }
        } /* end of thread interleaving loop */
        if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {
            break;
        } /* end of break condition check */
    } /* end of while(true) */
    candidate_thread_pool.clear();
    threads_to_exec.clear();
    progress.clear();
    /* USEING STATIC SCHEUDLING */
    engine = ChunkEngine(chunk_size, (2048 - 0));
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        candidate_thread_pool.emplace_back(tid_to_run);
    }
    while(true) {
        if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {
            while(!candidate_thread_pool.empty()) {
                tid_to_run = *(candidate_thread_pool.begin());
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] Assign chunk " << engine.getCurrentChunkRange() << endl;
#endif
                Chunk c = engine.getNextChunk(tid_to_run);
                Progress p("F_addr0", {c.first, 0}, c);
                progress[tid_to_run] = p;
                candidate_thread_pool.erase(candidate_thread_pool.begin());
                threads_to_exec.emplace_back(tid_to_run);
            } /* end of progress assignment */
        } /* end of chunk availability check */
        /* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
            if (!progress[tid_to_run].isInBound()) {
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] " << progress[tid_to_run].iteration[0] << " > " << progress[tid_to_run].chunk.second << endl;
#endif
                continue;
            }
#ifdef SIMULATOR_DEBUG
            cout << "[" << tid_to_run << "] Iterate " << progress[tid_to_run].ref << " at " << progress[tid_to_run].getIteration() << endl;
#endif
            if (progress[tid_to_run].ref == "F_addr0") {
                cnt++;
                int access = calAddrF_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_F_addr.find(access) != LAT_F_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_F_addr[access], 1.0);
                    }
                }
                LAT_F_addr[access] = cnt;
/* F_addr0 and C_addr0 are not in the same loop level */
/*   %k = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("C_addr0");
                continue;
            } /* end of check to F_addr0 */
            if (progress[tid_to_run].ref == "C_addr0") {
                cnt++;
                int access = calAddrC_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_C_addr.find(access) != LAT_C_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_C_addr[access], 1.0);
                    }
                }
                LAT_C_addr[access] = cnt;
                progress[tid_to_run].increment("D_addr0");
                continue;
            } /* end of check to C_addr0 */
            if (progress[tid_to_run].ref == "D_addr0") {
                cnt++;
                int access = calAddrD_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_D_addr.find(access) != LAT_D_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_D_addr[access], 1.0);
                    }
                }
                LAT_D_addr[access] = cnt;
                progress[tid_to_run].increment("F_addr1");
                continue;
            } /* end of check to D_addr0 */
            if (progress[tid_to_run].ref == "F_addr1") {
                cnt++;
                int access = calAddrF_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_F_addr.find(access) != LAT_F_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_F_addr[access], 1.0);
                    }
                }
                LAT_F_addr[access] = cnt;
                progress[tid_to_run].increment("F_addr2");
                continue;
            } /* end of check to F_addr1 */
            if (progress[tid_to_run].ref == "F_addr2") {
                cnt++;
                int access = calAddrF_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_F_addr.find(access) != LAT_F_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_F_addr[access], 1.0);
                    }
                }
                LAT_F_addr[access] = cnt;
/* F_addr2 is the last access node in a loop */
/*   %k = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[2] + 1 < 2048) {
                    progress[tid_to_run].iteration[2] += 1;
                    progress[tid_to_run].increment("C_addr0");
                    continue;
                } /* end of check to F_addr2 */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 2048) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("F_addr0");
                    continue;
                } /* end of check to F_addr2 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].iteration.pop_back();
/*   %j = alloca i32, align 4 has no depend parent */
                    progress[tid_to_run].iteration.emplace_back(0);
                    progress[tid_to_run].increment("F_addr0");
                    continue;
                } /* end of check to F_addr2 */
                if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {
                    candidate_thread_pool.emplace_back(tid_to_run);
                }
            }
        } /* end of thread interleaving loop */
        if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {
            break;
        } /* end of break condition check */
    } /* end of while(true) */
    candidate_thread_pool.clear();
    threads_to_exec.clear();
    progress.clear();
    /* USEING STATIC SCHEUDLING */
    engine = ChunkEngine(chunk_size, (2048 - 0));
    for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
        candidate_thread_pool.emplace_back(tid_to_run);
    }
    while(true) {
        if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {
            while(!candidate_thread_pool.empty()) {
                tid_to_run = *(candidate_thread_pool.begin());
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] Assign chunk " << engine.getCurrentChunkRange() << endl;
#endif
                Chunk c = engine.getNextChunk(tid_to_run);
                Progress p("G_addr0", {c.first, 0}, c);
                progress[tid_to_run] = p;
                candidate_thread_pool.erase(candidate_thread_pool.begin());
                threads_to_exec.emplace_back(tid_to_run);
            } /* end of progress assignment */
        } /* end of chunk availability check */
        /* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */
        for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {
            if (!progress[tid_to_run].isInBound()) {
#ifdef SIMULATOR_DEBUG
                cout << "[" << tid_to_run << "] " << progress[tid_to_run].iteration[0] << " > " << progress[tid_to_run].chunk.second << endl;
#endif
                continue;
            }
#ifdef SIMULATOR_DEBUG
            cout << "[" << tid_to_run << "] Iterate " << progress[tid_to_run].ref << " at " << progress[tid_to_run].getIteration() << endl;
#endif
            if (progress[tid_to_run].ref == "G_addr0") {
                cnt++;
                int access = calAddrG_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_G_addr.find(access) != LAT_G_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_G_addr[access], 1.0);
                    }
                }
                LAT_G_addr[access] = cnt;
/* G_addr0 and E_addr3 are not in the same loop level */
/*   %k = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("E_addr3");
                continue;
            } /* end of check to G_addr0 */
            if (progress[tid_to_run].ref == "E_addr3") {
                cnt++;
                int access = calAddrE_addr3(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_E_addr.find(access) != LAT_E_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_E_addr[access], 1.0);
                    }
                }
                LAT_E_addr[access] = cnt;
                progress[tid_to_run].increment("F_addr3");
                continue;
            } /* end of check to E_addr3 */
            if (progress[tid_to_run].ref == "F_addr3") {
                cnt++;
                int access = calAddrF_addr3(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_F_addr.find(access) != LAT_F_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_F_addr[access], 1.0);
                    }
                }
                LAT_F_addr[access] = cnt;
                progress[tid_to_run].increment("G_addr1");
                continue;
            } /* end of check to F_addr3 */
            if (progress[tid_to_run].ref == "G_addr1") {
                cnt++;
                int access = calAddrG_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_G_addr.find(access) != LAT_G_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_G_addr[access], 1.0);
                    }
                }
                LAT_G_addr[access] = cnt;
                progress[tid_to_run].increment("G_addr2");
                continue;
            } /* end of check to G_addr1 */
            if (progress[tid_to_run].ref == "G_addr2") {
                cnt++;
                int access = calAddrG_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1], progress[tid_to_run].iteration[2]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_G_addr.find(access) != LAT_G_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_G_addr[access], 1.0);
                    }
                }
                LAT_G_addr[access] = cnt;
/* G_addr2 is the last access node in a loop */
/*   %k = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[2] + 1 < 2048) {
                    progress[tid_to_run].iteration[2] += 1;
                    progress[tid_to_run].increment("E_addr3");
                    continue;
                } /* end of check to G_addr2 */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 2048) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("G_addr0");
                    continue;
                } /* end of check to G_addr2 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].iteration.pop_back();
/*   %j = alloca i32, align 4 has no depend parent */
                    progress[tid_to_run].iteration.emplace_back(0);
                    progress[tid_to_run].increment("G_addr0");
                    continue;
                } /* end of check to G_addr2 */
                if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {
                    candidate_thread_pool.emplace_back(tid_to_run);
                }
            }
        } /* end of thread interleaving loop */
        if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {
            break;
        } /* end of break condition check */
    } /* end of while(true) */
    candidate_thread_pool.clear();
    threads_to_exec.clear();
    progress.clear();
}
int main() {
#ifdef PAPI_TIMER
    PAPI_timer_init();
    PAPI_timer_start();
#endif
    interleaving();
    RTtoMR_AET();
#ifdef PAPI_TIMER
    PAPI_timer_end();
    PAPI_timer_print();
#endif
    rtDump();
    dumpMR();
    return 0;
}
 /* Analyze function: mm2 */ 
