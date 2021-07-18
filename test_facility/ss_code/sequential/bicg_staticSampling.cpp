
 /* Start to analysis array index
Array index info: Total number of references: 10
s.addr i
q.addr i
q.addr i
A.addr ((i * 8192) + j)
p.addr j
q.addr i
s.addr i
r.addr j
A.addr ((j * 8192) + i)
s.addr i
BC Array cost info: Total number of arrays: 5
A.addr 4
p.addr 2
q.addr 8
r.addr 2
s.addr 8
BC Average cost: 4.800000e+00

 Finish to analysis array index */ 

 /* Start to analyze argument
i32 %nx
i32 %ny
double* %A
double* %r
double* %s
double* %p
double* %q

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
Array r_addr 
Array p_addr 
Array A_addr 
Array q_addr 
Array s_addr 
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----array access s_addr i
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----array access q_addr i
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access q_addr i
------array access A_addr ((i * 8192) + j)
------array access p_addr j
------array access q_addr i
--i
--Loop Bound: (0, 8192)
--Loop inc: (i + 1)
--Loop predicate: <
----j
----Loop Bound: (0, 8192)
----Loop inc: (j + 1)
----Loop predicate: <
------array access s_addr i
------array access r_addr j
------array access A_addr ((j * 8192) + i)
------array access s_addr i

Finish analysis loops */ 
/* Start IV Dependence Analysis Finish to analyze IV Dependence *//* Start to analysis the access graphs_addr -> r_addr
r_addr -> A_addr
A_addr -> s_addr
s_addr -> s_addr
s_addr -> s_addr

q_addr -> q_addr
q_addr -> A_addr
A_addr -> p_addr
p_addr -> q_addr
q_addr -> q_addr
q_addr -> q_addr

s_addr -> s_addr

*/
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 819
----Sample number: 819
------Sample number: 671088
----Sample number: 819
------Sample number: 671088
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
std::unordered_map<uint64_t, uint64_t> LAT_s_addr;
std::unordered_map<uint64_t, uint64_t> LAT_q_addr;
std::unordered_map<uint64_t, uint64_t> LAT_A_addr;
std::unordered_map<uint64_t, uint64_t> LAT_p_addr;
std::unordered_map<uint64_t, uint64_t> LAT_r_addr;
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
/* s_addr i 0 */
int calAddrs_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* q_addr i 0 */
int calAddrq_addr0( int i) {
    int result = (i) * 8 / 64;
    return result;
}
/* q_addr i 1 */
int calAddrq_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* A_addr ((i * 8192) + j) 0 */
int calAddrA_addr0( int i, int j) {
    int result = (((i * 8192) + j)) * 8 / 64;
    return result;
}
/* p_addr j 0 */
int calAddrp_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* q_addr i 2 */
int calAddrq_addr2( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* s_addr i 1 */
int calAddrs_addr1( int i, int j) {
    int result = (i) * 8 / 64;
    return result;
}
/* r_addr j 0 */
int calAddrr_addr0( int i, int j) {
    int result = (j) * 8 / 64;
    return result;
}
/* A_addr ((j * 8192) + i) 1 */
int calAddrA_addr1( int i, int j) {
    int result = (((j * 8192) + i)) * 8 / 64;
    return result;
}
/* s_addr i 2 */
int calAddrs_addr2( int i, int j) {
    int result = (i) * 8 / 64;
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
    engine = ChunkEngine(chunk_size, (8192 - 0));
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
                Progress p("s_addr0", {c.first, }, c);
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
            if (progress[tid_to_run].ref == "s_addr0") {
                cnt++;
                int access = calAddrs_addr0(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_s_addr.find(access) != LAT_s_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_s_addr[access], 1.0);
                    }
                }
                LAT_s_addr[access] = cnt;
/* s_addr0 is the last access node in a loop */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].increment("s_addr0");
                    continue;
                } /* end of check to s_addr0 */
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
    engine = ChunkEngine(chunk_size, (8192 - 0));
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
                Progress p("q_addr0", {c.first, }, c);
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
            if (progress[tid_to_run].ref == "q_addr0") {
                cnt++;
                int access = calAddrq_addr0(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
/* q_addr0 and q_addr1 are not in the same loop level */
/*   %j = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("q_addr1");
                continue;
            } /* end of check to q_addr0 */
            if (progress[tid_to_run].ref == "q_addr1") {
                cnt++;
                int access = calAddrq_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
                progress[tid_to_run].increment("A_addr0");
                continue;
            } /* end of check to q_addr1 */
            if (progress[tid_to_run].ref == "A_addr0") {
                cnt++;
                int access = calAddrA_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_A_addr.find(access) != LAT_A_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_A_addr[access], 1.0);
                    }
                }
                LAT_A_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr0");
                continue;
            } /* end of check to A_addr0 */
            if (progress[tid_to_run].ref == "p_addr0") {
                cnt++;
                int access = calAddrp_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr2");
                continue;
            } /* end of check to p_addr0 */
            if (progress[tid_to_run].ref == "q_addr2") {
                cnt++;
                int access = calAddrq_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
/* q_addr2 is the last access node in a loop */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 8192) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("q_addr1");
                    continue;
                } /* end of check to q_addr2 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("q_addr0");
                    continue;
                } /* end of check to q_addr2 */
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
    engine = ChunkEngine(chunk_size, (8192 - 0));
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
                Progress p("s_addr1", {c.first, 0}, c);
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
            if (progress[tid_to_run].ref == "s_addr1") {
                cnt++;
                int access = calAddrs_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_s_addr.find(access) != LAT_s_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_s_addr[access], 1.0);
                    }
                }
                LAT_s_addr[access] = cnt;
                progress[tid_to_run].increment("r_addr0");
                continue;
            } /* end of check to s_addr1 */
            if (progress[tid_to_run].ref == "r_addr0") {
                cnt++;
                int access = calAddrr_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_r_addr.find(access) != LAT_r_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_r_addr[access], 1.0);
                    }
                }
                LAT_r_addr[access] = cnt;
                progress[tid_to_run].increment("A_addr1");
                continue;
            } /* end of check to r_addr0 */
            if (progress[tid_to_run].ref == "A_addr1") {
                cnt++;
                int access = calAddrA_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_A_addr.find(access) != LAT_A_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_A_addr[access], 1.0);
                    }
                }
                LAT_A_addr[access] = cnt;
                progress[tid_to_run].increment("s_addr2");
                continue;
            } /* end of check to A_addr1 */
            if (progress[tid_to_run].ref == "s_addr2") {
                cnt++;
                int access = calAddrs_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_s_addr.find(access) != LAT_s_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_s_addr[access], 1.0);
                    }
                }
                LAT_s_addr[access] = cnt;
/* s_addr2 is the last access node in a loop */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 8192) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("s_addr1");
                    continue;
                } /* end of check to s_addr2 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
/*   %j = alloca i32, align 4 has no depend parent */
                    progress[tid_to_run].iteration.emplace_back(0);
                    progress[tid_to_run].increment("s_addr1");
                    continue;
                } /* end of check to s_addr2 */
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
 /* Analyze function: bicg_cpu */ 
