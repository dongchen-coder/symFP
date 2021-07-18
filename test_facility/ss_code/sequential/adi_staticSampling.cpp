
 /* Start to analysis array index
Array index info: Total number of references: 34
v.addr ((0 + i) + 1)
p.addr (((i + 1) * 8192) + j)
u.addr (((i + 1) * 8192) + 0)
q.addr (((i + 1) * 8192) + 0)
u.addr (((j + 1) * 8192) + i)
p.addr (((i + 1) * 8192) + 0)
v.addr (0 + (i + 1))
q.addr (((i + 1) * 8192) + 0)
q.addr ((((i + 1) * 8192) + j) + 1)
v.addr (67100672 + (i + 1))
p.addr (((i + 1) * 8192) + (8190 - j))
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((8190 - j) * 8192) + (i + 1))
u.addr (((i + 1) * 8192) + 0)
p.addr (((i + 1) * 8192) + 0)
u.addr ((((j + 1) * 8192) + i) + 1)
u.addr ((((j + 1) * 8192) + i) + 2)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
q.addr ((((i + 1) * 8192) + j) + 1)
u.addr ((((i + 1) * 8192) + 8192) - 1)
p.addr (((((i + 1) * 8192) + 8192) - 2) - j)
u.addr (((((i + 1) * 8192) + 8192) - 1) - j)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
v.addr ((((8191 - j) * 8192) + i) + 1)
q.addr (((((i + 1) * 8192) + 8192) - 2) - j)
p.addr (((i + 1) * 8192) + j)
p.addr ((((i + 1) * 8192) + j) + 1)
v.addr (((i * 8192) + j) + 1)
v.addr ((((i + 1) * 8192) + j) + 1)
v.addr ((((i + 2) * 8192) + j) + 1)
q.addr (((i + 1) * 8192) + j)
p.addr (((i + 1) * 8192) + j)
u.addr (((((i + 1) * 8192) + 8192) - 2) - j)
BC Array cost info: Total number of arrays: 4
p.addr 24
q.addr 20
u.addr 19
v.addr 19
BC Average cost: 2.050000e+01

 Finish to analysis array index */ 

 /* Start to analyze argument
double* %p
double* %q
double* %v
double* %u

 Start to analysis argument */ 

 /* Start to analysis global variable 

 Finish to analysis global variable */ 

 /* Start analysis loops
Array q_addr 
Array u_addr 
Array p_addr 
Array v_addr 
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access v_addr ((0 + i) + 1)
----array access p_addr (((i + 1) * 8192) + 0)
----array access v_addr (0 + (i + 1))
----array access q_addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p_addr (((i + 1) * 8192) + j)
------array access p_addr ((((i + 1) * 8192) + j) + 1)
------array access u_addr (((j + 1) * 8192) + i)
------array access u_addr ((((j + 1) * 8192) + i) + 1)
------array access u_addr ((((j + 1) * 8192) + i) + 2)
------array access q_addr (((i + 1) * 8192) + j)
------array access p_addr (((i + 1) * 8192) + j)
------array access q_addr ((((i + 1) * 8192) + j) + 1)
----array access v_addr (67100672 + (i + 1))
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p_addr (((i + 1) * 8192) + (8190 - j))
------array access v_addr ((((8191 - j) * 8192) + i) + 1)
------array access q_addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access v_addr (((8190 - j) * 8192) + (i + 1))
--i
--Loop Bound: (0, 8190)
--Loop inc: (i + 1)
--Loop predicate: <
----array access u_addr (((i + 1) * 8192) + 0)
----array access p_addr (((i + 1) * 8192) + 0)
----array access u_addr (((i + 1) * 8192) + 0)
----array access q_addr (((i + 1) * 8192) + 0)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p_addr (((i + 1) * 8192) + j)
------array access p_addr ((((i + 1) * 8192) + j) + 1)
------array access v_addr (((i * 8192) + j) + 1)
------array access v_addr ((((i + 1) * 8192) + j) + 1)
------array access v_addr ((((i + 2) * 8192) + j) + 1)
------array access q_addr (((i + 1) * 8192) + j)
------array access p_addr (((i + 1) * 8192) + j)
------array access q_addr ((((i + 1) * 8192) + j) + 1)
----array access u_addr ((((i + 1) * 8192) + 8192) - 1)
----j
----Loop Bound: (0, 8190)
----Loop inc: (j + 1)
----Loop predicate: <
------array access p_addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u_addr (((((i + 1) * 8192) + 8192) - 1) - j)
------array access q_addr (((((i + 1) * 8192) + 8192) - 2) - j)
------array access u_addr (((((i + 1) * 8192) + 8192) - 2) - j)

Finish analysis loops */ 
/* Start IV Dependence Analysis Finish to analyze IV Dependence *//* Start to analysis the access graphu_addr -> p_addr
p_addr -> u_addr
u_addr -> q_addr
q_addr -> p_addr
p_addr -> p_addr
p_addr -> v_addr
v_addr -> v_addr
v_addr -> v_addr
v_addr -> q_addr
q_addr -> p_addr
p_addr -> q_addr
q_addr -> p_addr
q_addr -> u_addr
u_addr -> p_addr
p_addr -> u_addr
u_addr -> q_addr
q_addr -> u_addr
u_addr -> p_addr
u_addr -> u_addr

v_addr -> p_addr
p_addr -> v_addr
v_addr -> q_addr
q_addr -> p_addr
p_addr -> p_addr
p_addr -> u_addr
u_addr -> u_addr
u_addr -> u_addr
u_addr -> q_addr
q_addr -> p_addr
p_addr -> q_addr
q_addr -> p_addr
q_addr -> v_addr
v_addr -> p_addr
p_addr -> v_addr
v_addr -> q_addr
q_addr -> v_addr
v_addr -> p_addr
v_addr -> v_addr

*/
 /* Start to analysis the number of samples
calculating:
Dump tree:
----Sample number: 819
------Sample number: 670761
------Sample number: 670761
----Sample number: 819
------Sample number: 670761
------Sample number: 670761
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
std::unordered_map<uint64_t, uint64_t> LAT_v_addr;
std::unordered_map<uint64_t, uint64_t> LAT_p_addr;
std::unordered_map<uint64_t, uint64_t> LAT_u_addr;
std::unordered_map<uint64_t, uint64_t> LAT_q_addr;
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
/* v_addr ((0 + i) + 1) 0 */
int calAddrv_addr0( int i) {
    int result = (((0 + i) + 1)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + 0) 0 */
int calAddrp_addr0( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* v_addr (0 + (i + 1)) 1 */
int calAddrv_addr1( int i) {
    int result = ((0 + (i + 1))) * 8 / 64;
    return result;
}
/* q_addr (((i + 1) * 8192) + 0) 0 */
int calAddrq_addr0( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + j) 1 */
int calAddrp_addr1( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* p_addr ((((i + 1) * 8192) + j) + 1) 2 */
int calAddrp_addr2( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* u_addr (((j + 1) * 8192) + i) 0 */
int calAddru_addr0( int i, int j) {
    int result = ((((j + 1) * 8192) + i)) * 8 / 64;
    return result;
}
/* u_addr ((((j + 1) * 8192) + i) + 1) 1 */
int calAddru_addr1( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* u_addr ((((j + 1) * 8192) + i) + 2) 2 */
int calAddru_addr2( int i, int j) {
    int result = (((((j + 1) * 8192) + i) + 2)) * 8 / 64;
    return result;
}
/* q_addr (((i + 1) * 8192) + j) 1 */
int calAddrq_addr1( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + j) 3 */
int calAddrp_addr3( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* q_addr ((((i + 1) * 8192) + j) + 1) 2 */
int calAddrq_addr2( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* v_addr (67100672 + (i + 1)) 2 */
int calAddrv_addr2( int i) {
    int result = ((67100672 + (i + 1))) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + (8190 - j)) 4 */
int calAddrp_addr4( int i, int j) {
    int result = ((((i + 1) * 8192) + (8190 - j))) * 8 / 64;
    return result;
}
/* v_addr ((((8191 - j) * 8192) + i) + 1) 3 */
int calAddrv_addr3( int i, int j) {
    int result = (((((8191 - j) * 8192) + i) + 1)) * 8 / 64;
    return result;
}
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 3 */
int calAddrq_addr3( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* v_addr (((8190 - j) * 8192) + (i + 1)) 4 */
int calAddrv_addr4( int i, int j) {
    int result = ((((8190 - j) * 8192) + (i + 1))) * 8 / 64;
    return result;
}
/* u_addr (((i + 1) * 8192) + 0) 3 */
int calAddru_addr3( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + 0) 5 */
int calAddrp_addr5( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* u_addr (((i + 1) * 8192) + 0) 4 */
int calAddru_addr4( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* q_addr (((i + 1) * 8192) + 0) 4 */
int calAddrq_addr4( int i) {
    int result = ((((i + 1) * 8192) + 0)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + j) 6 */
int calAddrp_addr6( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* p_addr ((((i + 1) * 8192) + j) + 1) 7 */
int calAddrp_addr7( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* v_addr (((i * 8192) + j) + 1) 5 */
int calAddrv_addr5( int i, int j) {
    int result = ((((i * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* v_addr ((((i + 1) * 8192) + j) + 1) 6 */
int calAddrv_addr6( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* v_addr ((((i + 2) * 8192) + j) + 1) 7 */
int calAddrv_addr7( int i, int j) {
    int result = (((((i + 2) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* q_addr (((i + 1) * 8192) + j) 5 */
int calAddrq_addr5( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* p_addr (((i + 1) * 8192) + j) 8 */
int calAddrp_addr8( int i, int j) {
    int result = ((((i + 1) * 8192) + j)) * 8 / 64;
    return result;
}
/* q_addr ((((i + 1) * 8192) + j) + 1) 6 */
int calAddrq_addr6( int i, int j) {
    int result = (((((i + 1) * 8192) + j) + 1)) * 8 / 64;
    return result;
}
/* u_addr ((((i + 1) * 8192) + 8192) - 1) 5 */
int calAddru_addr5( int i) {
    int result = (((((i + 1) * 8192) + 8192) - 1)) * 8 / 64;
    return result;
}
/* p_addr (((((i + 1) * 8192) + 8192) - 2) - j) 9 */
int calAddrp_addr9( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* u_addr (((((i + 1) * 8192) + 8192) - 1) - j) 6 */
int calAddru_addr6( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 1) - j)) * 8 / 64;
    return result;
}
/* q_addr (((((i + 1) * 8192) + 8192) - 2) - j) 7 */
int calAddrq_addr7( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
    return result;
}
/* u_addr (((((i + 1) * 8192) + 8192) - 2) - j) 7 */
int calAddru_addr7( int i, int j) {
    int result = ((((((i + 1) * 8192) + 8192) - 2) - j)) * 8 / 64;
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
    engine = ChunkEngine(chunk_size, (8190 - 0));
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
                Progress p("v_addr0", {c.first, }, c);
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
            if (progress[tid_to_run].ref == "v_addr0") {
                cnt++;
                int access = calAddrv_addr0(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr0");
                continue;
            } /* end of check to v_addr0 */
            if (progress[tid_to_run].ref == "p_addr0") {
                cnt++;
                int access = calAddrp_addr0(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr1");
                continue;
            } /* end of check to p_addr0 */
            if (progress[tid_to_run].ref == "v_addr1") {
                cnt++;
                int access = calAddrv_addr1(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr0");
                continue;
            } /* end of check to v_addr1 */
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
/* q_addr0 and p_addr1 are not in the same loop level */
/*   %j = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("p_addr1");
                continue;
            } /* end of check to q_addr0 */
            if (progress[tid_to_run].ref == "p_addr1") {
                cnt++;
                int access = calAddrp_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr2");
                continue;
            } /* end of check to p_addr1 */
            if (progress[tid_to_run].ref == "p_addr2") {
                cnt++;
                int access = calAddrp_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr0");
                continue;
            } /* end of check to p_addr2 */
            if (progress[tid_to_run].ref == "u_addr0") {
                cnt++;
                int access = calAddru_addr0(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr1");
                continue;
            } /* end of check to u_addr0 */
            if (progress[tid_to_run].ref == "u_addr1") {
                cnt++;
                int access = calAddru_addr1(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr2");
                continue;
            } /* end of check to u_addr1 */
            if (progress[tid_to_run].ref == "u_addr2") {
                cnt++;
                int access = calAddru_addr2(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr1");
                continue;
            } /* end of check to u_addr2 */
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
                progress[tid_to_run].increment("p_addr3");
                continue;
            } /* end of check to q_addr1 */
            if (progress[tid_to_run].ref == "p_addr3") {
                cnt++;
                int access = calAddrp_addr3(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
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
            } /* end of check to p_addr3 */
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
                if (progress[tid_to_run].iteration[1] + 1 < 8190) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("p_addr1");
                    continue;
                } /* end of check to q_addr2 */
                progress[tid_to_run].iteration.pop_back();
                progress[tid_to_run].increment("v_addr2");
                continue;
            } /* end of check to q_addr2 */
            if (progress[tid_to_run].ref == "v_addr2") {
                cnt++;
                int access = calAddrv_addr2(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
/* v_addr2 and p_addr4 are not in the same loop level */
/*   %j = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("p_addr4");
                continue;
            } /* end of check to v_addr2 */
            if (progress[tid_to_run].ref == "p_addr4") {
                cnt++;
                int access = calAddrp_addr4(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr3");
                continue;
            } /* end of check to p_addr4 */
            if (progress[tid_to_run].ref == "v_addr3") {
                cnt++;
                int access = calAddrv_addr3(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr3");
                continue;
            } /* end of check to v_addr3 */
            if (progress[tid_to_run].ref == "q_addr3") {
                cnt++;
                int access = calAddrq_addr3(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr4");
                continue;
            } /* end of check to q_addr3 */
            if (progress[tid_to_run].ref == "v_addr4") {
                cnt++;
                int access = calAddrv_addr4(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
/* v_addr4 is the last access node in a loop */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 8190) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("p_addr4");
                    continue;
                } /* end of check to v_addr4 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("v_addr0");
                    continue;
                } /* end of check to v_addr4 */
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
    engine = ChunkEngine(chunk_size, (8190 - 0));
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
                Progress p("u_addr3", {c.first, }, c);
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
            if (progress[tid_to_run].ref == "u_addr3") {
                cnt++;
                int access = calAddru_addr3(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr5");
                continue;
            } /* end of check to u_addr3 */
            if (progress[tid_to_run].ref == "p_addr5") {
                cnt++;
                int access = calAddrp_addr5(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr4");
                continue;
            } /* end of check to p_addr5 */
            if (progress[tid_to_run].ref == "u_addr4") {
                cnt++;
                int access = calAddru_addr4(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr4");
                continue;
            } /* end of check to u_addr4 */
            if (progress[tid_to_run].ref == "q_addr4") {
                cnt++;
                int access = calAddrq_addr4(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
/* q_addr4 and p_addr6 are not in the same loop level */
/*   %j = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("p_addr6");
                continue;
            } /* end of check to q_addr4 */
            if (progress[tid_to_run].ref == "p_addr6") {
                cnt++;
                int access = calAddrp_addr6(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr7");
                continue;
            } /* end of check to p_addr6 */
            if (progress[tid_to_run].ref == "p_addr7") {
                cnt++;
                int access = calAddrp_addr7(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr5");
                continue;
            } /* end of check to p_addr7 */
            if (progress[tid_to_run].ref == "v_addr5") {
                cnt++;
                int access = calAddrv_addr5(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr6");
                continue;
            } /* end of check to v_addr5 */
            if (progress[tid_to_run].ref == "v_addr6") {
                cnt++;
                int access = calAddrv_addr6(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("v_addr7");
                continue;
            } /* end of check to v_addr6 */
            if (progress[tid_to_run].ref == "v_addr7") {
                cnt++;
                int access = calAddrv_addr7(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_v_addr.find(access) != LAT_v_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_v_addr[access], 1.0);
                    }
                }
                LAT_v_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr5");
                continue;
            } /* end of check to v_addr7 */
            if (progress[tid_to_run].ref == "q_addr5") {
                cnt++;
                int access = calAddrq_addr5(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
                progress[tid_to_run].increment("p_addr8");
                continue;
            } /* end of check to q_addr5 */
            if (progress[tid_to_run].ref == "p_addr8") {
                cnt++;
                int access = calAddrp_addr8(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr6");
                continue;
            } /* end of check to p_addr8 */
            if (progress[tid_to_run].ref == "q_addr6") {
                cnt++;
                int access = calAddrq_addr6(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
/* q_addr6 is the last access node in a loop */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 8190) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("p_addr6");
                    continue;
                } /* end of check to q_addr6 */
                progress[tid_to_run].iteration.pop_back();
                progress[tid_to_run].increment("u_addr5");
                continue;
            } /* end of check to q_addr6 */
            if (progress[tid_to_run].ref == "u_addr5") {
                cnt++;
                int access = calAddru_addr5(progress[tid_to_run].iteration[0]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
/* u_addr5 and p_addr9 are not in the same loop level */
/*   %j = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration.emplace_back(0);
                progress[tid_to_run].increment("p_addr9");
                continue;
            } /* end of check to u_addr5 */
            if (progress[tid_to_run].ref == "p_addr9") {
                cnt++;
                int access = calAddrp_addr9(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_p_addr.find(access) != LAT_p_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_p_addr[access], 1.0);
                    }
                }
                LAT_p_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr6");
                continue;
            } /* end of check to p_addr9 */
            if (progress[tid_to_run].ref == "u_addr6") {
                cnt++;
                int access = calAddru_addr6(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
                progress[tid_to_run].increment("q_addr7");
                continue;
            } /* end of check to u_addr6 */
            if (progress[tid_to_run].ref == "q_addr7") {
                cnt++;
                int access = calAddrq_addr7(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_q_addr.find(access) != LAT_q_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_q_addr[access], 1.0);
                    }
                }
                LAT_q_addr[access] = cnt;
                progress[tid_to_run].increment("u_addr7");
                continue;
            } /* end of check to q_addr7 */
            if (progress[tid_to_run].ref == "u_addr7") {
                cnt++;
                int access = calAddru_addr7(progress[tid_to_run].iteration[0], progress[tid_to_run].iteration[1]);
                /* choose a number between 1 and 100 */ 
                int enable = rand() % 100 + 1;
                if (enable <= 10) {
                    if (LAT_u_addr.find(access) != LAT_u_addr.end()) {
                        rtHistoCal(RT, cnt - LAT_u_addr[access], 1.0);
                    }
                }
                LAT_u_addr[access] = cnt;
/* u_addr7 is the last access node in a loop */
/*   %j = alloca i32, align 4 has no depend parent */
                if (progress[tid_to_run].iteration[1] + 1 < 8190) {
                    progress[tid_to_run].iteration[1] += 1;
                    progress[tid_to_run].increment("p_addr9");
                    continue;
                } /* end of check to u_addr7 */
/*   %i = alloca i32, align 4 has no depend parent */
                progress[tid_to_run].iteration[0] += 1;
                if (progress[tid_to_run].isInBound()) {
                    progress[tid_to_run].iteration.pop_back();
                    progress[tid_to_run].increment("u_addr3");
                    continue;
                } /* end of check to u_addr7 */
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
 /* Analyze function: adi */ 
