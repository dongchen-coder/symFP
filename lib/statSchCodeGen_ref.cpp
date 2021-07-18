#include "statSchCodeGen_ref.hpp"

namespace statSchCodeGen_ref {
	char StaticSchedulingSamplerCodeGen::ID = 0;
	static RegisterPass<StaticSchedulingSamplerCodeGen> X("statSchCodeGen_ref", "parallel random interleaving static sampling code generating pass (both reference and iteration)", false, false);

	StaticSchedulingSamplerCodeGen::StaticSchedulingSamplerCodeGen() : FunctionPass(ID) {}

	void StaticSchedulingSamplerCodeGen::numberRefToSameArray(LoopRefTNode *LoopRefTree) {
		
		if (LoopRefTree->AA != NULL) {
			if (refToSameArrayCnt.find(arrayName[LoopRefTree->AA]) == refToSameArrayCnt.end()) {
				refNumber[LoopRefTree->AA] = 0;
				refToSameArrayCnt[arrayName[LoopRefTree->AA]] = 1;
			} else {
				refNumber[LoopRefTree->AA] = refToSameArrayCnt[arrayName[LoopRefTree->AA]];
				refToSameArrayCnt[arrayName[LoopRefTree->AA]] += 1;
			}
		}
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				numberRefToSameArray(*it);
			}
		}
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::numberLoops(LoopRefTNode *LoopRefTree) {
		
		if (LoopRefTree->L != NULL) {
			if (loopNumber.find(LoopRefTree->L) == loopNumber.end()) {
				loopNumber[LoopRefTree->L] = loopNumber.size();
			}
		}
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				numberLoops(*it);
			}
		}
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::initIndvName(LoopRefTNode *LoopRefTree) {
		
		if (LoopRefTree == NULL) {
			return;
		}
		
		if (LoopRefTree->L != NULL) {
			for (std::vector<Value *>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
				indvName[*it] = (*it)->getName().str();
				if (std::find(indvName[*it].begin(), indvName[*it].end(), '.') != indvName[*it].end()) {
					indvName[*it].erase(std::find(indvName[*it].begin(), indvName[*it].end(), '.'));
				}
			}
		}
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode *>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				initIndvName(*it);
			}
		}
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::initArrayName() {
		
		for (std::map<Instruction*, std::string>::iterator it = arrayName.begin(), eit = arrayName.end(); it != eit; ++it) {
			it->second.replace(std::find(it->second.begin(), it->second.end(), '.'), std::find(it->second.begin(), it->second.end(), '.') +1, 1, '_');
		}
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::addrCalFuncGen(LoopRefTNode* LoopRefTree, std::vector<string> indvs) {
		
		if (LoopRefTree->L != NULL) {
			for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
				indvs.push_back(indvName[(*it)]);
			}
		}
		
		if (LoopRefTree->AA != NULL) {
			errs() << "/* " + arrayName[LoopRefTree->AA] + " " + arrayExpression[LoopRefTree->AA] + " " + std::to_string(refNumber[LoopRefTree->AA])  + " */\n";
			errs() << "int calAddr" + arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
			string arguments = "";
			for (std::vector<string>::iterator it = indvs.begin(), eit = indvs.end(); it != eit; ++it) {
				arguments += "int " + (*it) + ", ";
			}
			if (arguments.size() != 0) {
				arguments.pop_back();
				arguments.pop_back();
			}
			errs() << arguments;
			errs() << ") {\n";
#if defined (CLS) && defined (DS)
			errs() << "    int result = (" + arrayExpression[LoopRefTree->AA] + ") * " + std::to_string(DS) + " / " + std::to_string(CLS) + ";\n";
#else
			errs() << "    int result = " + arrayExpression[LoopRefTree->AA] + ";\n";
#endif
			errs() << "    return result;\n";
			errs() << "}\n";
			return;
		}
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				addrCalFuncGen(*it, indvs);
			}
		}
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::addrCalFuncGenTop(LoopRefTNode *LoopRefTree) {
		
		std::vector<string> indvs;
		addrCalFuncGen(LoopRefTree, indvs);
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::headerGen() {
		
		// errs() << "#define DEBUG\n";
		
		errs() << "#include <map>\n";
		errs() << "#include <set>\n";
		errs() << "#include <vector>\n";
		errs() << "#include <tuple>\n";
		errs() << "#include <unordered_map>\n";
		errs() << "#include <algorithm>\n";
		errs() << "#include <stdlib.h>\n";
		errs() << "#include <iostream>\n";
		errs() << "#include <iomanip>\n";
		errs() << "#include <cmath>\n";
		errs() << "#include <time.h>\n";
		errs() << "#include <cassert>\n";
#ifdef PARALLEL_CXX_THREAD
		errs() << "#include <thread>\n";
		errs() << "#include <mutex>\n";
#endif
		errs() << "#ifdef PAPI_TIMER\n";
		errs() << "#  include \"papi_timer.h\"\n";
		errs() << "#endif\n";
#ifdef PARALLEL
		errs() << "#ifndef THREAD_NUM\n";
		errs() << "#    define THREAD_NUM   4\n";
		errs() << "#endif\n";
		errs() << "#ifndef BIN_SIZE\n";
		errs() << "#    define BIN_SIZE   4\n";
		errs() << "#endif\n";
		errs() << "#ifndef CHUNK_SIZE\n";
		errs() << "#    define CHUNK_SIZE   4\n";
		errs() << "#endif\n";
#endif

		errs() << "using namespace std;\n";
		errs() << "typedef pair<int, int> Chunk;\n";
#ifdef PARALLEL_CXX_THREAD
		errs() << "std::mutex mtx;\n";
#endif
#ifdef REFERENCE_GROUP
		errs() << "std::map<string, map<uint64_t, double>> refRT;\n";
#endif
		errs() << "class ChunkEngine {\n";
		errs() << "    int lb = 0;\n";
		errs() << "    int ub = 0;\n";
		errs() << "    int chunk_size = 0;\n";
		errs() << "    int trip = 0;\n";
		errs() << "    int avail_chunk = 0;\n";
		errs() << "public:\n";
		errs() << "    ChunkEngine() {} \n";
		errs() << "    ChunkEngine(int chunk_size, int trip) {\n";
		errs() << "        assert(chunk_size <= trip);\n";
		errs() << "        this->chunk_size = chunk_size;\n";
		errs() << "        this->trip = trip;\n";
		errs() << "        this->avail_chunk = (trip / chunk_size + (trip % chunk_size));\n";
		errs() << "        this->lb = 0;\n";
		errs() << "        this->ub = chunk_size - 1;\n";
		errs() << "    }\n";
		errs() << "    string getCurrentChunkRange() {\n";
		errs() << "        return \"[\" + to_string(this->lb) + \", \" + to_string(this->ub) + \"]\";\n";
		errs() << "    }\n";
		errs() << "    bool hasNextChunk() {\n";
		errs() << "        return this->avail_chunk > 0;\n";
		errs() << "    }\n";
		errs() << "    Chunk getNextChunk(int tid) {\n";
		errs() << "        // assign the current lb, ub to thread tid and update the next chunk\n";
		errs() << "        Chunk curr = make_pair(this->lb, this->ub);\n";
		errs() << "        this->lb = this->ub + 1;\n";
		errs() << "        this->ub = (this->lb + chunk_size - 1) <= this->trip ? (this->lb + chunk_size - 1) : this->trip;\n";
		errs() << "        this->avail_chunk -= 1;\n";
		errs() << "        return curr;\n";
		errs() << "    }\n";
		errs() << "};\n";
		errs() << "class Progress {\n";
		errs() << "public:\n";
		errs() << "    string ref;\n";
		errs() << "    Chunk chunk;\n";
		errs() << "    vector<int> iteration;\n";
		errs() << "    Progress() { }\n";
		errs() << "    Progress(string ref, vector<int> iteration, Chunk c) {\n";
		errs() << "        this->ref = ref;\n";
		errs() << "        this->iteration = iteration;\n";
		errs() << "        this->chunk = c;\n";
		errs() << "    }\n";
		errs() << "    string getIteration() {\n";
		errs() << "        string ret = \"(\";\n";
		errs() << "        for (int i = 0; i < this->iteration.size(); i++) {\n";
		errs() << "            ret += to_string(this->iteration[i]);\n";
		errs() << "            if (i != this->iteration.size() - 1)\n";
		errs() << "                ret += \",\";\n";
		errs() << "            }\n";
		errs() << "        ret += \")\";\n";
		errs() << "        return ret;\n";
		errs() << "    }\n";
		errs() << "    string getReference() {\n";
		errs() << "        return this->ref;\n";
		errs() << "    }\n";
		errs() << "    void increment(string ref, vector<int> iteration) {\n";
		errs() << "        this->ref = ref;\n";
		errs() << "        this->iteration = iteration;\n";
		errs() << "    }\n";
		errs() << "    void increment(string ref) {\n";
		errs() << "        this->ref = ref;\n";
		errs() << "    }\n";
		errs() << "    bool isInBound() {\n";
		errs() << "        assert(this->iteration[0] >= chunk.first);\n";
		errs() << "        return this->iteration[0] <= chunk.second;\n";
		errs() << "    }\n";
		errs() << "};\n";
		//errs() << "std::map<uint64_t, tuple<uint64_t, int>> LAT;\n";
		vector<string> visit;
		for (std::map<Instruction*, int>::iterator it = refNumber.begin(), eit = refNumber.end(); it != eit; ++it) {
			if (find(visit.begin(), visit.end(), arrayName[it->first]) != visit.end()) { continue; }
			visit.push_back(arrayName[it->first]);
			errs() << "std::unordered_map<uint64_t, uint64_t> LAT_" + arrayName[it->first] + ";\n";
		}
		
		errs() << "std::map<uint64_t, double> RT;\n";
		errs() << "std::map<uint64_t, double> MR;\n";
		return;
	}

#ifdef DumpRTMR
	void StaticSchedulingSamplerCodeGen::rtHistoGen() {
		errs() << "void rtHistoCal( map<uint64_t, double> &rth, int rt, int val ) {\n";
		errs() << "    if ( val <= 0) {\n";
		errs() << "        return;\n";
		errs() << "    }\n";
#ifdef PARALLEL_OMP
		errs() << "    #pragma omp critical\n";
		errs() << "    {";
		errs() << "        if (rth.find(rt) == rth.end()) { \n";
		errs() << "            rth[rt] = val;\n";
		errs() << "        } else {\n";
		errs() << "            rth[rt] += val;\n";
		errs() << "        }\n";
		errs() << "    }";
#elif defined(PARALLEL_CXX_THREAD)
//		errs() << "    std::unique_lock<std::mutex> lck (mtx,std::defer_lock);\n";
//		errs() << "    lck.lock();\n";
		errs() << "    if (rth.find(rt) == rth.end()) { \n";
		errs() << "        rth[rt] = val;\n";
		errs() << "    } else {\n";
		errs() << "        rth[rt] += val;\n";
		errs() << "    }\n";
//		errs() << "    lck.unlock();\n";
#else
		errs() << "    if (rth.find(rt) == rth.end()) { \n";
		errs() << "        rth[rt] = val;\n";
		errs() << "    } else {\n";
		errs() << "        rth[rt] += val;\n";
		errs() << "    }\n";
#endif
		errs() << "    return;\n";
		errs() << "}\n";
#ifdef REFERENCE_GROUP   
		errs() << "\n";    
		errs() << "void refRTHistoCal( int rt, int val, string ref ) {\n";
		errs() << "    if ( val <= 0) {\n;";
		errs() << "        return;\n";
		errs() << "    }\n";
		errs() << "    if (refRT.find(ref) == refRT.end()) { \n";
		errs() << "        refRT[ref][rt] = val;\n";
		errs() << "    }\n";
		errs() << "    else {\n";
		errs() << "        if (refRT[ref].find(rt) == refRT[ref].end()) { \n";
		errs() << "            refRT[ref][rt] = val;\n";
		errs() << "        } else {\n";
		errs() << "            refRT[ref][rt] += val;\n";
		errs() << "        }\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
#endif
	}
#endif

	/* Generate the function to calculate the bins */
	/* Generate the function to calculate the bins */
	void StaticSchedulingSamplerCodeGen::subBlkRTGen() {
		string space = "    ";
		errs() << "void subBlkRT(map<uint64_t, double> &rth, int rt) {\n";
#ifdef PARALLEL_CXX_THREAD
		errs() << space + "unique_lock< mutex> lck (mtx, defer_lock);\n";
		errs() << space + "lck.lock();\n";
#endif
		errs() << space + "int msb = 0;\n";
		errs() << space + "int tmp_rt = rt;\n";
		errs() << space + "while(tmp_rt != 0) {\n";
		errs() << space + "    tmp_rt = tmp_rt / 2;\n";
		errs() << space + "    ++msb;\n";
		errs() << space + "}\n";
		errs() << space + "if (msb >= BIN_SIZE) {\n";
		errs() << space + "    int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;\n";
		errs() << space + "    for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {\n";
		errs() << space + "        if (rt < b) {\n";
		errs() << space + "            rtHistoCal(rth, b - diff, 1);\n";
		errs() << space + "            break;\n";
		errs() << space + "        }\n";
		errs() << space + "    }\n";
		errs() << space + "}\n";
		errs() << space + "else {\n";
		errs() << space + "    rtHistoCal(rth, pow(2, msb-1), 1);\n";
		errs() << space + "}\n";
#ifdef PARALLEL_CXX_THREAD
		errs() << space + "lck.unlock();\n";
#endif
		errs() << space + "return;\n";
		errs() << "}\n";
#ifdef REFERENCE_GROUP  
		errs() << "\n";
		errs() << "void refSubBlkRT(int rt, string ref) {\n";
		errs() << space + "int msb = 0;\n";
		errs() << space + "int tmp_rt = rt;\n";
		errs() << space + "while(tmp_rt != 0) {\n";
		errs() << space + "    tmp_rt = tmp_rt / 2;\n";
		errs() << space + "    ++msb;\n";
		errs() << space + "}\n";
		errs() << space + "if (msb >= BIN_SIZE) {\n";
		errs() << space + "    int diff = (pow(2, msb) - pow(2, msb-1)) / BIN_SIZE;\n";
		errs() << space + "    for (int b = pow(2, msb-1); b <= pow(2, msb); b+=diff) {\n";
		errs() << space + "        if (rt < b) {\n";
		errs() << space + "            refRTHistoCal(b - diff, 1, ref);\n";
		errs() << space + "            break;\n";
		errs() << space + "        }\n";
		errs() << space + "    }\n";
		errs() << space + "}\n";
		errs() << space + "else {\n";
		errs() << space + "    refRTHistoCal(pow(2, msb-1), 1, ref);\n";
		errs() << space + "}\n";
		errs() << space + "return;\n";
		errs() << "}\n";
#endif
		return;
	}
	
	
#ifdef DumpRTMR
	void StaticSchedulingSamplerCodeGen::rtDumpGen() {
		
		errs() << "void rtDump() {\n";
		errs() << "    cout << \"Start to dump reuse time histogram\\n\";\n";
		errs() << "    cout << fixed << setprecision(3);\n";
		errs() << "    for (map<uint64_t, double>::iterator it = RT.begin(), eit = RT.end(); it != eit; ++it) {\n";
		errs() << "        cout << it->first << \", \" << it->second << \"\\n\";\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
#ifdef REFERENCE_GROUP
		errs() << "\n";
		errs() << "void refRTDump() {\n";
		errs() << "    for (map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
		errs() << "        cout << \"Start to dump reuse time histogram for \" << it->first << \"\\n\";\n";
		errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(), eiit = it->second.end(); iit != eiit; ++iit) {\n";
		errs() << "            cout << it->first << \",\" << iit->first << \",\" << iit->second << endl;\n";
		errs() << "        }\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
#endif  
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::rtToMRGen() {
		
		errs() << "void RTtoMR_AET() {\n";
		
		string space = "    ";
		errs() << space + "std::map<uint64_t, double> P;\n";
		errs() << space + "double total_num_RT = 0;\n";
		errs() << space + "uint64_t max_RT = 0;\n";
		
		errs() << space + "for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
		errs() << space + "    total_num_RT += it->second;\n";
		errs() << space + "    if (max_RT < it->first) {\n";
		errs() << space + "        max_RT = it->first;\n";
		errs() << space + "    }\n";
		errs() << space + "}\n";
		
		errs() << space + "double accumulate_num_RT = 0;\n";
		errs() << space + "for (std::map<uint64_t, double>::reverse_iterator it = RT.rbegin(), eit = RT.rend(); it != eit; ++it) {\n";
		errs() << space + "    P[it->first] = accumulate_num_RT / total_num_RT;\n";
		errs() << space + "    accumulate_num_RT += it->second;\n";
		errs() << space + "}\n";
		
		errs() << space + "P[0] = 1;\n";
		
		errs() << space + "double sum_P = 0;\n";
		errs() << space + "uint64_t t = 0;\n";
		errs() << space + "uint64_t prev_t = 0;\n";
//        errs() << space + "for (uint64_t c = 0; c <= max_RT; c++) {\n";
// 128M cache, 32B cache line
//        errs() << space + "for (uint64_t c = 0; c <= max_RT && c <= 1677216; c++) {\n";
// 20MB
		errs() << space + "for (uint64_t c = 0; c <= max_RT && c <= 327680; c++) {\n";
		errs() << space + "    while (sum_P < c && t <= max_RT) {\n";
		errs() << space + "        if (P.find(t) != P.end()) {\n";
		errs() << space + "            sum_P += P[t];\n";
		errs() << space + "            prev_t = t;\n";
		errs() << space + "        } else {\n";
		errs() << space + "            sum_P += P[prev_t];\n";
		errs() << space + "        }\n";
		errs() << space + "        t++;\n";
		errs() << space + "    }\n";
		
		errs() << space + "    MR[c] = P[prev_t];\n";
		errs() << space + "}\n";
		
		errs() << space + "return;\n";
		
		errs() << "}\n";
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::mrDumpGen() {
		
		errs() << "void dumpMR() {\n";
		
		errs() << "    cout << \"miss ratio\" << endl;\n";
		
		errs() << "    std::map<uint64_t, double>::iterator it1 = MR.begin();\n";
		errs() << "    std::map<uint64_t, double>::iterator it2 = MR.begin();\n";
		
		errs() << "    while(it1 != MR.end()) {\n";
		errs() << "        while(1) {\n";
		errs() << "            std::map<uint64_t, double>::iterator it3 = it2;\n";
		errs() << "            ++it3;\n";
		errs() << "            if (it3 == MR.end()) {\n";
		errs() << "                break;\n";
		errs() << "            }\n";
		errs() << "            if (it1->second - it3->second < 0.00001) {\n";
		errs() << "                ++it2;\n";
		errs() << "            } else {\n";
		errs() << "                break;\n";
		errs() << "            }\n";
		errs() << "        }\n";
		errs() << "        cout << it1->first << \", \" << it1->second << endl;\n";
		errs() << "        if (it1 != it2) {\n";
		errs() << "            cout << it2->first << \", \" << it2->second << endl;\n";
		errs() << "        }\n";
		errs() << "        it1 = ++it2;\n";
		errs() << "        it2 = it1;\n";
		errs() << "    }\n";
		
		errs() << "    return;\n";
		errs() << "}\n";
		
		return;
	}
#endif

#if defined(PARALLEL) && defined(REFERENCE_GROUP)
	void StaticSchedulingSamplerCodeGen::rtMergeGen() {
		errs() << "/* Merge the refRT to RT */\n";
		errs() << "void rtMerge() {\n";
		errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
		errs() << "        for (map<uint64_t, double>::iterator iit = it->second.begin(); iit != it->second.end(); ++iit) {\n"; 
		errs() << "            rtHistoCal(iit->first, iit->second);\n";
		errs() << "        }\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
		return;
	}

	void StaticSchedulingSamplerCodeGen::GaussianDistrGen() {
		errs() << "/* Distribute the per-reference RTHisto based on Gasussian */\n";
		errs() << "void gaussian_distr() {\n";
		errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
		errs() << "        map<uint64_t, double> tmp = it->second;\n";
		errs() << "        for (map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {\n"; 
		errs() << "            uint64_t mu = iit->first;\n";
		errs() << "            uint64_t Y = iit->second;\n";
		errs() << "            double sigma = 3 * (double)mu / 7;\n";
		errs() << "            uint64_t start_b = log2(mu) <= log2(THREAD_NUM) ? 0 : (log2(mu) - log2(THREAD_NUM));\n";
		errs() << "            uint64_t end_b = log2(mu) + log2(THREAD_NUM);\n";
		errs() << "            /* Clear the original refRT[ref]. Make sure that this will be done only once */\n";
		errs() << "            if (iit == tmp.begin()) {\n";
		errs() << "                refRT[it->first].clear();\n";
		errs() << "            }\n";
		errs() << "            for(int b = start_b; b <= end_b; b++) {\n";
		errs() << "                if (b >= 0) {\n";
		errs() << "                    double val = 0.0;\n";
		errs() << "                    val = (1 /  (sqrt(2 * M_PI) * sigma))* exp( -1 * pow((pow(2.0,(double)b) - mu), 2.0) / (2 * pow(sigma, 2.0)));\n";
		errs() << "                    cout << pow(2.0,(double)b) << \", \" << (val * Y) << endl;\n";
		errs() << "                    refRTHistoCal(pow(2.0,(double)b), (val * Y), it->first);\n";
		errs() << "                }\n";
		errs() << "            }\n";
		errs() << "        }\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
		return;
		
	}

	void StaticSchedulingSamplerCodeGen::UniformDistrGen() {
		errs() << "/* Distribute the per-reference RTHisto Uniformly. Equally splite the RT to a range of bins */\n";
		errs() << "void uniform_distr() {\n";
		errs() << "    for(map<string, map<uint64_t, double>>::iterator it = refRT.begin(); it != refRT.end(); ++it) {\n";
		errs() << "        map<uint64_t, double> tmp = it->second;\n";
		errs() << "        for (map<uint64_t, double>::iterator iit = tmp.begin(); iit != tmp.end(); ++iit) {\n";
		errs() << "            uint64_t mu = iit->first;\n";
		errs() << "            uint64_t start_b = log2(mu) <= log2(THREAD_NUM) ? 0 : (log2(mu) - log2(THREAD_NUM));\n";
		errs() << "            uint64_t end_b = (log2(mu) + log2(THREAD_NUM));\n";
		errs() << "            uint64_t split_val = iit->second / (end_b - start_b + 1);\n";
		errs() << "            /* Clear the original refRT[ref]. Make sure that this will be done only once */\n";
		errs() << "            if (iit == tmp.begin()) {\n";
		errs() << "                refRT[it->first].clear();\n";
		errs() << "            }\n";
		errs() << "            for(int b = start_b; b <= end_b; b++) {\n";
		errs() << "                if (b >= 0) {\n";
		errs() << "                    refRTHistoCal(pow(2.0, b), split_val, it->first);\n";
		errs() << "                }\n";
		errs() << "            }\n";
		errs() << "        }\n";
		errs() << "    }\n";
		errs() << "    return;\n";
		errs() << "}\n";
		return;
	}
#endif

	string StaticSchedulingSamplerCodeGen::getBound(Value *bound) {
		
		if (isa<Instruction>(bound)) {
			
			Instruction *inst = cast<Instruction>(bound);
			
			switch (inst->getOpcode()) {
				case Instruction::Add:
					return "(" + getBound(inst->getOperand(0)) + " + " + getBound(inst->getOperand(1)) + ")";
					break;
				case Instruction::Sub:
					return "(" + getBound(inst->getOperand(0)) + " - " + getBound(inst->getOperand(1)) + ")";;
					break;
				case Instruction::Mul:
					return "(" + getBound(inst->getOperand(0)) + " * " + getBound(inst->getOperand(1)) + ")";;
					break;
				case Instruction::FDiv:
				case Instruction::SDiv:
				case Instruction::UDiv:
					return "(" + getBound(inst->getOperand(0)) + " / " + getBound(inst->getOperand(1)) + ")";;
					break;
				case Instruction::Load:
					return inst->getOperand(0)->getName().str();
					break;
				case Instruction::Alloca:
					return inst->getName().str();
				default:
					break;
			}
		}
		else if (isa<ConstantInt>(bound)) {
			return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
		}
		return "";
	}
	
	string StaticSchedulingSamplerCodeGen::getBound_Start(Value *bound) {
		if (isa<Instruction>(bound)) {
			
			Instruction *inst = cast<Instruction>(bound);
			
			switch (inst->getOpcode()) {
				case Instruction::Add:
					return "(" + getBound_Start(inst->getOperand(0)) + " + " + getBound_Start(inst->getOperand(1)) + ")";
					break;
				case Instruction::Sub:
					return "(" + getBound_Start(inst->getOperand(0)) + " - " + getBound_Start(inst->getOperand(1)) + ")";;
					break;
				case Instruction::Mul:
					return "(" + getBound_Start(inst->getOperand(0)) + " * " + getBound_Start(inst->getOperand(1)) + ")";;
					break;
				case Instruction::FDiv:
				case Instruction::SDiv:
				case Instruction::UDiv:
					return "(" + getBound_Start(inst->getOperand(0)) + " / " + getBound_Start(inst->getOperand(1)) + ")";;
					break;
				case Instruction::Load:
					return inst->getOperand(0)->getName().str() + "_Start";
					break;
				case Instruction::Alloca:
					return inst->getName().str() + "_Start";
				default:
					break;
			}
		}
		else if (isa<ConstantInt>(bound)) {
			return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
		}
		return "";
	}
	
	std::vector<LoopRefTNode*> StaticSchedulingSamplerCodeGen::findLoops(
		LoopRefTNode *LoopRefTree,
		std::vector<LoopRefTNode*> ref
	) {
		
		/*
		if (LoopRefTree->L != NULL) {
			loops.push_back(LoopRefTree);
		}
		
		if (LoopRefTree->AA != NULL && LoopRefTree == ref) {
			return loops;
		}
		
		
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				std::vector<LoopRefTNode*> loopTmp = findLoops(*it, loops, enableOPT);
				if (loopTmp.size() != 0) {
					loopRes = loopTmp;
				}
			}
		}
		 */
		std::vector<LoopRefTNode*> loopRes;
		return loopRes;
	}

	DepParentInfo StaticSchedulingSamplerCodeGen::getParent(Value * v)
	{
		int depType = 0; // 1 means LB dep, -1 means UB dep
		DepNode * parent = getAnalysis<ivdepAnalysis::IVDependenceAnalysis>().getParentOfValue(v, &depType);
		return std::make_pair(parent, depType);
	}

#if 0
	/* Generate the loop iteration vector incremental func */
	void StaticSchedulingSamplerCodeGen::LoopIterIncGen(AccGraph *G,
															AccGraphEdge *E,
															LoopRefTNode* IDomLoop,
															string space)
	{
		string upper_bound = getBound(IDomLoop->LIS->LB->front().second);
		DepParentInfo pinfo = getParent(IDomLoop->LIS->LB->front().second);
		if (pinfo.first && pinfo.second < 0)
			upper_bound = "progress[tid_to_run].iteration[" + to_string(pinfo.second-1) + "]";
		string progressRepr = "progress[tid_to_run].iteration[" + to_string(IDomLoop->LoopLevel-1) + "]";
		errs() << space << progressRepr;
		errs() << space << "if (" << progressRepr << " + 1 " << getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().predicateToString(IDomLoop->LIS->PREDICATE->front()) << " " << upper_bound << ") {\n";
		errs() << space << "    " << progressRepr << " += 1;\n";
		errs() << space << "    continue; /* go back to the interleaving */ \n";
		errs() << space << "}\n";
		return;
	}
#endif

	/* Generate the loop iteration vector update func */
	void StaticSchedulingSamplerCodeGen::LoopIterUpdateGen(AccGraph *G,
															AccGraphEdge *E,
															LoopRefTNode* IDomLoop,
															string space,
															bool EnableIterationInc,
															bool EnableIterationUpdate)
	{
		// increment part
		if (EnableIterationInc) {
			string upper_bound = getBound(IDomLoop->LIS->LB->front().second);
			DepParentInfo pinfo = getParent(IDomLoop->LIS->IDV->front());
			if (pinfo.first && pinfo.second < 0)
				upper_bound = "progress[tid_to_run].iteration[" + to_string(pinfo.first->level-1) + "]";
			// this is the outermost loop
			string progressRepr = "progress[tid_to_run].iteration[" + to_string(IDomLoop->LoopLevel-1) + "]";
			if (IDomLoop->LoopLevel-1 == 0) {
				upper_bound = "progress[tid_to_run].chunk.second + 1";
				errs() << space << progressRepr << " += 1;\n";
				errs() << space << "if (progress[tid_to_run].isInBound()) {\n";
			} else {
				errs() << space << "if (" << progressRepr << " + 1 " << getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().predicateToString(IDomLoop->LIS->PREDICATE->front()) << " " <<  upper_bound << ") {\n";
				errs() << space << "    " << progressRepr << " += 1;\n";
			}
		}
		if (EnableIterationUpdate) {
			// update part
			// pop iteration vector L1.size() times
			// push LB of loop in L2 in iteration vector
			vector<LoopRefTNode *> L1, L2;
			vector<LoopRefTNode *>::iterator loopIter;
			G->QueryEngine.getPath(E->getSourceNode(), IDomLoop, L1);
			G->QueryEngine.getPath(E->getSinkNode(), IDomLoop, L2);
			for (loopIter = L1.begin(); loopIter != L1.end(); ++loopIter) {
				errs() << space << "    " << "progress[tid_to_run].iteration.pop_back();\n";
			}
			for (loopIter = L2.begin(); loopIter != L2.end(); ++loopIter) {
				string lower_bound = getBound((*loopIter)->LIS->LB->front().first);
				DepParentInfo pinfo = getParent((*loopIter)->LIS->IDV->front());
				if (pinfo.first && pinfo.second > 0) {
					lower_bound = "progress[tid_to_run].iteration[" + to_string(pinfo.first->level - 1) + "]";
				}
				errs() << space << "    " << "progress[tid_to_run].iteration.emplace_back(" << lower_bound << ");\n";
			}
		}
//		if (EnableIterationInc) {
//			errs() << space << "    continue; /* go back to the interleaving */ \n";
//			errs() << space << "}\n";
//		}
		
		return;
	}

	void StaticSchedulingSamplerCodeGen::reuseUpdateGen(LoopRefTNode *access, int currentAccessLoopNestCnt, string accessName, string space)
	{
		errs() << space << "cnt++;\n";
		errs() << space << "int access = calAddr" << accessName << "(";
		for(int i = 0; i < currentAccessLoopNestCnt; i++) {
			errs() << "progress[tid_to_run].iteration[" << i << "]";
			if (i != currentAccessLoopNestCnt-1) {
				errs() << ", ";
			}
		}
		errs() << ");\n";
		errs() << space << "/* choose a number between 1 and 100 */ \n";
#ifdef PARALLEL_CXX_THREAD
		errs() << space << "std::unique_lock<std::mutex> lck (mtx,std::defer_lock);\n";
		errs() << space << "lck.lock();\n";
#endif
		errs() << space << "int enable = rand() % 100 + 1;\n";
		errs() << space << "if (enable <= " << to_string((int)(SamplingRate * 100)) << ") {\n";
		errs() << space << "    if (LAT_" << arrayName[access->AA] << ".find(access) != LAT_" << arrayName[access->AA] << ".end()) {\n";
		errs() << space << "        rtHistoCal(RT, cnt - LAT_" << arrayName[access->AA] << "[access], 1.0);\n";
		errs() << space << "    }\n";
		errs() << space << "}\n";
#ifdef PARALLEL_CXX_THREAD
		errs() << space << "lck.unlock();\n";
#endif
		errs() << space << "LAT_" << arrayName[access->AA] << "[access] = cnt;\n";
	}
	

	/* Generating Reuse Search */
	bool StaticSchedulingSamplerCodeGen::refRTSearchGen(
		LoopRefTNode *root,
		LoopRefTNode *LoopRefTree,
		vector<LoopRefTNode*> currentLoops,
		string space
		) {
		/**
		 Generate the preparation code at the very beginning of the outermost loop
		 */
		if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
			currentLoops.emplace_back(LoopRefTree);
			AccGraph * G = getAnalysis<AccGraphAnalysis::AccessGraphAnalysis>().AccGraphLoopMap[LoopRefTree];
			/* If there is no array reference in current LoopNest, we do not need to generate
			 any code.
			 */
			vector<LoopRefTNode *> NodesInGraph = G->getAllGraphNodes();
			vector<LoopRefTNode *>::iterator GraphIter = NodesInGraph.begin();
			LoopRefTNode * currentNode = NULL, * nextNode = NULL;
			bool reachTheLastAccessNode = false;
			while (true) {
				currentNode = (*GraphIter);
				vector<LoopRefTNode *> pathToNode = G->getPathToNode(currentNode, NULL);
				// generate the preparion work at the very beginning
				// this will be generate if and only if currently we are in the root
				// of the access graph
				if (GraphIter == NodesInGraph.begin()) {
					errs() << space << "engine = ChunkEngine(chunk_size, (" << getBound(LoopRefTree->LIS->LB->front().second) << " - " << getBound(LoopRefTree->LIS->LB->front().first) << "));\n";
					errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
					errs() << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
					errs() << space << "}\n";
					errs() << space << "while(true) {\n";
#if defined(RANDOM_INTERLEAVING)
					errs() << "chunk_assignment:\n";
#endif
					errs() << space << space << "if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {\n";
					errs() << space << space << space << "while(!candidate_thread_pool.empty()) {\n";
					errs() << space << space << space << space << "tid_to_run = *(candidate_thread_pool.begin());\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] Assign chunk \" << engine.getCurrentChunkRange() << endl;\n";
					errs() << "#endif\n";
					errs() << space << space << space << space << "Chunk c = engine.getNextChunk(tid_to_run);\n";
					errs() << space << space << space << space << "Progress p(\"" + arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", {";
			
					vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
					// traverse the loop from outer to inner
					while (pathRIter != pathToNode.rend()) {
						if (pathRIter == pathToNode.rbegin()) { //
							errs() << "c.first, ";
						} else {
							DepParentInfo pinfo = getParent((*pathRIter)->LIS->IDV->front());
							string lower_bound = getBound((*pathRIter)->LIS->LB->front().first);
							if (pinfo.first && pinfo.second > 0) { // has LB dependence
								lower_bound = "progress[tid_to_run].iteration[" + to_string(pinfo.first->level - 1) + "]";
							}
							errs() << lower_bound;
							if ( pathRIter != pathToNode.rend()-1)
								errs() << ", ";
						}
						pathRIter++;
					}
					errs() << "}, c);\n";
					errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
					errs() << space << space << space << space << "candidate_thread_pool.erase(tid_to_run);\n";
					errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
					errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#endif
					errs() << space << space << space << "} /* end of progress assignment */\n";
					errs() << space << space << "} /* end of chunk availability check */\n";
					/* Generate the Interleaving code for the current outermost loop */
#ifdef UNIFORM_INTERLEAVING
					errs() << space << space << "/* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */\n";
//					errs() << space << space << "/* UNIFORM THREAD INTERLEAVING */\n";
					errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//					errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; });\n";
//					errs() << space << space << "for (auto tid_to_run : threads_to_exec) {\n";
#elif defined(RANDOM_INTERLEAVING)
					errs() << space << space << "/* RANDOM THREAD INTERLEAVING */\n";
					errs() << space << space << "while ( !threads_to_exec.empty()) {\n";
					errs() << space << space << space << "tid_to_run = threads_to_exec[rand() % threads_to_exec.size()];\n";
#endif
					errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] \" << progress[tid_to_run].iteration[0] << \" > \" << progress[tid_to_run].chunk.second << endl;\n";
					errs() << "#endif\n";
					errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
					errs() << space << space << space << space << "continue;\n";
					errs() << space << space << space << "}\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
					errs() << "#endif\n";
				}
				string currAccessName = arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
				errs() << space << space << space << "if (progress[tid_to_run].ref == \"" << currAccessName << "\") {\n";
				reuseUpdateGen(currentNode, (int)pathToNode.size(), currAccessName, space + space + space + space);
				// move to the next access
				GraphIter++;
				if (GraphIter == NodesInGraph.end()) {
					// this is the last node in the access graph
					// the next access node to exec is the header node
					GraphIter = NodesInGraph.begin();
					reachTheLastAccessNode = true;
				}
				nextNode = (*GraphIter);
				nextAccessName = arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
				// check whether we need to insert the increment funciton
				//
				// currentNode stores the current access to be examined, nextNode stores
				// the next access to be examined.
				//
				// If this access is the last access of a loop, then we need to generate the
				// loop increment function. In this case, we will query all edges that sourced
				// by the currentNode.
				//
				// Then there are three cases:
				// CASE 1:
				// Edge (currentNode -> XXX) is carried by a loop (backedge) and this loop is not
				// the immediate loop dominator of currentNode
				//
				// We need to compute loops betwen currentNode and the carried loop, L1;
				// loops betwen XXX and the carried loop, L2.
				// Then the iteration vector will pop L1.size times, and push the LB of all
				// loops in L2.
				// We also need to generate the loop IV increment code for IV of carried loop
				//
				// CASE 2:
				// Edge (currentNode -> nextNode)
				//
				// This has two cases:
				// a) currentNode is the last access node
				// b) others
				//
				// For a): We need to firstly compute the immediate dominator that dominates both
				// currentNode and nextNode. Then we do the same as a): compute loops between
				// currentNode and the common immediate dominator, L1; compute loops between nextNode
				// and the common immediate dominator, L2. Pop the iteration vector L1.size times,
				// push the LB of all loops in L2.
				//
				// For b): The same as a) but there is one more step, increment the LIV of the outermost
				// loop.
				//
				// CASE 3:
				// Edge (currentNode -> XXX) is carried by a loop (backedge) and this loop is the
				// immediate loop dominator of currentNode.
				//
				// increment the LIV of the carried loop
				//
				if (G->QueryEngine.isLastAccessNodeInALoop(currentNode) || !G->QueryEngine.areTwoAccessInSameLevel(currentNode, nextNode)) {
					if (G->QueryEngine.isLastAccessNodeInALoop(currentNode)) {
						errs() << "/* " << currAccessName << " is the last access node in a loop */\n";
					} else if (!G->QueryEngine.areTwoAccessInSameLevel(currentNode, nextNode)) {
						errs() << "/* " << currAccessName << " and " << nextAccessName << " are not in the same loop level */\n";
					}
					vector<AccGraphEdge *> edges = G->getEdgesWithSource(currentNode);
					LoopRefTNode * sink = NULL, * carryLoop = NULL;
					for (auto edge : edges) {
						carryLoop = edge->getCarryLoop();
						sink = edge->getSinkNode();
						nextAccessName = arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
						if (carryLoop) { // CASE 1 or CASE 3
							// CASE 3
							if (G->QueryEngine.getImmediateLoopDominator(currentNode) == carryLoop) {
								// increment the LIV only
								LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, false);
								errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
								errs() << space << space << space << space << space << "continue;\n";
								errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
							} else { // CASE 1
								// pop iteration vector L1.size() times
								// push LB of loop in L2 in iteration vector
								// and increment the LIV
								LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, true);
								errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
								errs() << space << space << space << space << space << "continue;\n";
								errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
							}
						} else if (edge->hasSource(currentNode) && edge->hasSink(nextNode)) { // CASE 2
							vector<LoopRefTNode *> dominator;
							G->QueryEngine.getImmediateCommonLoopDominator(currentNode, sink, reachTheLastAccessNode, dominator);
							// pop iteration vector L1.size() times
							// push LB of loop in L2 in iteration vector
							// increment the LIV if and only currentNode is the last access node
							LoopIterUpdateGen(G, edge, dominator.front(), space + space + space, reachTheLastAccessNode, true);
							errs() << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
							errs() << space << space << space << space << "continue;\n";
							errs() << space << space << space << "} /* end of check to " << currAccessName << " */\n";
						}
					}
				} else {
					errs() << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
					errs() << space << space << space << space << "continue;\n";
					errs() << space << space << space << "} /* end of check to " << currAccessName << " */\n";
				}
				if (reachTheLastAccessNode)
					break;
			}
//			errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
//			errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
			errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
			errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
			errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
//			errs() << space << space << space << space << "}\n";
			errs() << space << space << space << "}\n";
			errs() << space << space << "} /* end of thread interleaving loop */\n";
//			errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
			errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
			errs() << space << space << space << "break;\n";
			errs() << space << space << "} /* end of break condition check */\n";
			errs() << space << "} /* end of while(true) */\n";
		}
		return true;
	}


	void StaticSchedulingSamplerCodeGen::refRTBodyGen(LoopRefTNode *LoopRefTree) {
		
		string space = "    ";
		errs() << space << "uint64_t cnt = 0;\n";
		errs() << space << "int tid_to_run, chunk_size;\n";
		errs() << space << "set<int> candidate_thread_pool;\n";
		errs() << space << "vector<int> threads_to_exec;\n";
		errs() << space << "auto randgen = []() { return rand() % THREAD_NUM; };\n";
#if defined(RANDOM_INTERLEAVING)
		errs() << space << "srand((unsigned)time(NULL));\n";
#endif
		errs() << space << "ChunkEngine engine;\n";
		errs() << space << "unordered_map<int, Progress> progress;\n";
//		errs() << space << "uint64_t access;\n";
		errs() << "#ifdef CHUNK_SIZE\n";
		errs() << space << space << "chunk_size = CHUNK_SIZE;\n";
		errs() << "#else\n";
		errs() << space << space << "chunk_size = 1;\n";
		errs() << "#endif\n";
		if (LoopRefTree->L) {
			refRTSearchGen(LoopRefTree, LoopRefTree, vector<LoopRefTNode *>(), "    ");
		} else {
			/* We remove the dummy entry node and run the RTSearchCode for each outermost loop */
			for (auto subLoop : *(LoopRefTree->next)) {
				refRTSearchGen(subLoop, subLoop, vector<LoopRefTNode *>(), "    ");
				/* When change to another outermost loop, we should remove all metadata */
				errs() << space << "candidate_thread_pool.clear();\n";
				errs() << space << "threads_to_exec.clear();\n";
				errs() << space << "progress.clear();\n";
			}
		}
		return;
	}

	void StaticSchedulingSamplerCodeGen::perArrayRTBodyGen(std::string array, LoopRefTNode *LoopRefTree)
	{
		string space = "    ";
		errs() << space << "uint64_t cnt = 0;\n";
		errs() << space << "int tid_to_run, chunk_size;\n";
		errs() << space << "set<int> candidate_thread_pool;\n";
		errs() << space << "vector<int> threads_to_exec;\n";
		errs() << space << "auto randgen = []() { return rand() % THREAD_NUM; };\n";
#if defined(RANDOM_INTERLEAVING)
		errs() << space << "srand((unsigned)time(NULL));\n";
#endif
		errs() << space << "ChunkEngine engine;\n";
		errs() << space << "unordered_map<int, Progress> progress;\n";
//		errs() << space << "uint64_t access;\n";
		errs() << "#ifdef CHUNK_SIZE\n";
		errs() << space << space << "chunk_size = CHUNK_SIZE;\n";
		errs() << "#else\n";
		errs() << space << space << "chunk_size = 1;\n";
		errs() << "#endif\n";
		if (LoopRefTree->L) {
			refRTSearchGen(LoopRefTree, LoopRefTree, vector<LoopRefTNode *>(), "    ");
		} else {
			// We remove the dummy entry node and run the RTSearchCode for each outermost loop
			for (auto subLoop : *(LoopRefTree->next)) {
				// filter out those outermost loop that does not have the given array access
				if (!getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().DoesLoopContainsArray(subLoop, array))
					continue;
				perArrayRTSearchGen(array, subLoop, subLoop, vector<LoopRefTNode *>(), "    ");
				// When change to another outermost loop, we should remove all metadata
				errs() << space << "candidate_thread_pool.clear();\n";
				errs() << space << "threads_to_exec.clear();\n";
				errs() << space << "progress.clear();\n";
			}
		}
		return;
		
	}

	bool StaticSchedulingSamplerCodeGen::perArrayRTSearchGen(std::string array, LoopRefTNode *root, LoopRefTNode *LoopRefTree, std::vector<LoopRefTNode *> currentLoops, std::string space)
	{
		/**
		 Generate the preparation code at the very beginning of the outermost loop
		 */
		if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
			currentLoops.emplace_back(LoopRefTree);
			AccGraph * G = getAnalysis<AccGraphAnalysis::AccessGraphAnalysis>().AccGraphLoopMap[LoopRefTree];
			/* If there is no array reference in current LoopNest, we do not need to generate
			 any code.
			 */
			vector<LoopRefTNode *> NodesInGraph = G->getAllGraphNodes();
			vector<LoopRefTNode *>::iterator GraphIter = NodesInGraph.begin();
			vector<LoopRefTNode *>::iterator TargetArrayIter = NodesInGraph.begin();
			LoopRefTNode * currentNode = NULL, * nextNode = NULL;
			LoopRefTNode * currentTargetArrayNode = NULL, * nextTargetArrayNode = NULL;
			bool reachTheLastAccessNode = false;
			while (true) {
				currentNode = (*GraphIter);
				vector<LoopRefTNode *> pathToNode = G->getPathToNode(currentNode, NULL);
				// generate the preparion work at the very beginning
				// this will be generate if and only if currently we are in the root
				// of the access graph
				if (GraphIter == NodesInGraph.begin()) {
					errs() << space << "engine = ChunkEngine(chunk_size, (" << getBound(LoopRefTree->LIS->LB->front().second) << " - " << getBound(LoopRefTree->LIS->LB->front().first) << "));\n";
					errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
					errs() << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
					errs() << space << "}\n";
					errs() << space << "while(true) {\n";
#if defined(RANDOM_INTERLEAVING)
					errs() << "chunk_assignment:\n";
#endif
					errs() << space << space << "if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {\n";
					errs() << space << space << space << "while(!candidate_thread_pool.empty()) {\n";
					errs() << space << space << space << space << "tid_to_run = *(candidate_thread_pool.begin());\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] Assign chunk \" << engine.getCurrentChunkRange() << endl;\n";
					errs() << "#endif\n";
					errs() << space << space << space << space << "Chunk c = engine.getNextChunk(tid_to_run);\n";
					errs() << space << space << space << space << "Progress p(\"" + arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", {";
			
					vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
					// traverse the loop from outer to inner
					while (pathRIter != pathToNode.rend()) {
						if (pathRIter == pathToNode.rbegin()) { //
							errs() << "c.first, ";
						} else {
							DepParentInfo pinfo = getParent((*pathRIter)->LIS->IDV->front());
							string lower_bound = getBound((*pathRIter)->LIS->LB->front().first);
							if (pinfo.first && pinfo.second > 0) { // has LB dependence
								lower_bound = "progress[tid_to_run].iteration[" + to_string(pinfo.first->level - 1) + "]";
							}
							errs() << lower_bound;
							if ( pathRIter != pathToNode.rend()-1)
								errs() << ", ";
						}
						pathRIter++;
					}
					errs() << "}, c);\n";
					errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
					errs() << space << space << space << space << "candidate_thread_pool.erase(tid_to_run);\n";
					errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
					errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#endif
					errs() << space << space << space << "} /* end of progress assignment */\n";
					errs() << space << space << "} /* end of chunk availability check */\n";
					/* Generate the Interleaving code for the current outermost loop */
#ifdef UNIFORM_INTERLEAVING
					errs() << space << space << "/* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */\n";
//					errs() << space << space << "/* UNIFORM THREAD INTERLEAVING */\n";
					errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//					errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; });\n";
//					errs() << space << space << "for (auto tid_to_run : threads_to_exec) {\n";
#elif defined(RANDOM_INTERLEAVING)
					errs() << space << space << "/* RANDOM THREAD INTERLEAVING */\n";
					errs() << space << space << "while ( !threads_to_exec.empty()) {\n";
					errs() << space << space << space << "tid_to_run = threads_to_exec[rand() % threads_to_exec.size()];\n";
#endif
					errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] \" << progress[tid_to_run].iteration[0] << \" > \" << progress[tid_to_run].chunk.second << endl;\n";
					errs() << "#endif\n";
					errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
					errs() << space << space << space << space << "continue;\n";
					errs() << space << space << space << "}\n";
					errs() << "#ifdef SIMULATOR_DEBUG\n";
					errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
					errs() << "#endif\n";
				}
				string currAccessName = arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
				errs() << space << space << space << "if (progress[tid_to_run].ref == \"" << currAccessName << "\") {\n";
				if (arrayName[currentNode->AA] == array) {
					reuseUpdateGen(currentNode, (int)pathToNode.size(), currAccessName, space + space + space + space);
				} else {
					errs() << space << space << space << space << "cnt++; // skip " << currAccessName << "\n";
				}
				// move to the next access
				GraphIter++;
				if (GraphIter == NodesInGraph.end()) {
					// this is the last node in the access graph
					// the next access node to exec is the header node
					GraphIter = NodesInGraph.begin();
					reachTheLastAccessNode = true;
				}
				nextNode = (*GraphIter);
				// The purpose of this process if to find the next access node that
				// access the target array, skip those array accesses not relevant
				vector<LoopRefTNode *>::iterator tmpGraphIter = GraphIter;
				nextAccessName = arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
				
				// check whether we need to insert the increment funciton
				//
				// currentNode stores the current access to be examined, nextNode stores
				// the next access to be examined.
				//
				// If this access is the last access of a loop, then we need to generate the
				// loop increment function. In this case, we will query all edges that sourced
				// by the currentNode.
				//
				// Then there are three cases:
				// CASE 1:
				// Edge (currentNode -> XXX) is carried by a loop (backedge) and this loop is not
				// the immediate loop dominator of currentNode
				//
				// We need to compute loops betwen currentNode and the carried loop, L1;
				// loops betwen XXX and the carried loop, L2.
				// Then the iteration vector will pop L1.size times, and push the LB of all
				// loops in L2.
				// We also need to generate the loop IV increment code for IV of carried loop
				//
				// CASE 2:
				// Edge (currentNode -> nextNode)
				//
				// This has two cases:
				// a) currentNode is the last access node
				// b) others
				//
				// For a): We need to firstly compute the immediate dominator that dominates both
				// currentNode and nextNode. Then we do the same as a): compute loops between
				// currentNode and the common immediate dominator, L1; compute loops between nextNode
				// and the common immediate dominator, L2. Pop the iteration vector L1.size times,
				// push the LB of all loops in L2.
				//
				// For b): The same as a) but there is one more step, increment the LIV of the outermost
				// loop.
				//
				// CASE 3:
				// Edge (currentNode -> XXX) is carried by a loop (backedge) and this loop is the
				// immediate loop dominator of currentNode.
				//
				// increment the LIV of the carried loop
				//
				if (G->QueryEngine.isLastAccessNodeInALoop(currentNode) || !G->QueryEngine.areTwoAccessInSameLevel(currentNode, nextNode)) {
					if (G->QueryEngine.isLastAccessNodeInALoop(currentNode)) {
						errs() << "/* " << currAccessName << " is the last access node in a loop */\n";
					} else if (!G->QueryEngine.areTwoAccessInSameLevel(currentNode, nextNode)) {
						errs() << "/* " << currAccessName << " and " << nextAccessName << " are not in the same loop level */\n";
					}
					vector<AccGraphEdge *> edges = G->getEdgesWithSource(currentNode);
					LoopRefTNode * sink = NULL, * carryLoop = NULL;
					for (auto edge : edges) {
						carryLoop = edge->getCarryLoop();
						sink = edge->getSinkNode();
						nextAccessName = arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
						if (carryLoop) { // CASE 1 or CASE 3
							// CASE 3
							if (G->QueryEngine.getImmediateLoopDominator(currentNode) == carryLoop) {
								// increment the LIV only
								LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, false);
								errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
								errs() << space << space << space << space << space << "continue;\n";
								errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
							} else { // CASE 1
								// pop iteration vector L1.size() times
								// push LB of loop in L2 in iteration vector
								// and increment the LIV
								LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, true);
								errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
								errs() << space << space << space << space << space << "continue;\n";
								errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
							}
						} else if (edge->hasSource(currentNode) && edge->hasSink(nextNode)) { // CASE 2
							vector<LoopRefTNode *> dominator;
							G->QueryEngine.getImmediateCommonLoopDominator(currentNode, sink, reachTheLastAccessNode, dominator);
							// pop iteration vector L1.size() times
							// push LB of loop in L2 in iteration vector
							// increment the LIV if and only currentNode is the last access node
							LoopIterUpdateGen(G, edge, dominator.front(), space + space + space, reachTheLastAccessNode, true);
							errs() << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
							errs() << space << space << space << space << "continue;\n";
							errs() << space << space << space << "} /* end of check to " << currAccessName << " */\n";
						}
					}
				} else {
					errs() << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
					errs() << space << space << space << space << "continue;\n";
					errs() << space << space << space << "} /* end of check to " << currAccessName << " */\n";
				}
				if (reachTheLastAccessNode)
					break;
			}
//			errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
//			errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
			errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
			errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
			errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
//			errs() << space << space << space << space << "}\n";
			errs() << space << space << space << "}\n";
			errs() << space << space << "} /* end of thread interleaving loop */\n";
//			errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
			errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
			errs() << space << space << space << "break;\n";
			errs() << space << space << "} /* end of break condition check */\n";
			errs() << space << "} /* end of while(true) */\n";
		}
		return true;
	}
	
	void StaticSchedulingSamplerCodeGen::refRTGen(LoopRefTNode *LoopRefTree) {
//#ifdef PERREFERENCE
		std::unordered_set<std::string>::iterator arrayIter = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet.begin();
		for (; arrayIter != getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet.end(); ++arrayIter) {
			errs() << "void interleaving_" << *arrayIter << "() {\n";
			perArrayRTBodyGen(*arrayIter, LoopRefTree);
			errs() << "}\n";
		}
//#else
//		errs() << "void interleaving() {\n";
//		refRTBodyGen(LoopRefTree);
//		errs() << "}\n";
//#endif
		
		return;
	}
	
	void StaticSchedulingSamplerCodeGen::mainGen() {
	
		errs() << "int main() {\n";
		
		string space = "    ";
		// add timer
		errs() << "#ifdef PAPI_TIMER\n";
		errs() << space << "PAPI_timer_init();\n";
		errs() << space << "PAPI_timer_start();\n";
		errs() << "#endif\n";
		
		std::unordered_set<std::string>::iterator arrayIter = ArrayNameSet.begin();
		for (; arrayIter != ArrayNameSet.end(); ++arrayIter) {
#if defined(PARALLEL_CXX_THREAD)
			errs() << space + "thread t_"+ *arrayIter + "(";
			errs() << "interleaving_" + *arrayIter + ");\n";
#else
			errs() << space << "interleaving_" << *arrayIter << "();\n";
#endif
		}
#if defined(PARALLEL_CXX_THREAD)
		arrayIter = ArrayNameSet.begin();
		for (; arrayIter != ArrayNameSet.end(); ++arrayIter) {
			errs() << space + "t_"+ *arrayIter + ".join();\n";
		}
#endif
//		errs() << space << "interleaving();\n";

		errs() << "    RTtoMR_AET();\n";
		errs() << "#ifdef PAPI_TIMER\n";
		errs() << "    PAPI_timer_end();\n";
		errs() << "    PAPI_timer_print();\n";
		errs() << "#endif\n";

#ifdef DumpRTMR
#ifdef REFERENCE_GROUP
		// errs() << "    gaussian_distr();\n";
		// errs() << "    uniform_distr();\n";
		errs() << "    refRTDump();\n";  
		errs() << "    rtMerge();\n";
#endif
		errs() << "    rtDump();\n";
		errs() << "    dumpMR();\n";
#elif defined(DumpRefLease)
		errs() << "    RL_main(0);\n";
#endif
		
		errs() << "    return 0;\n";
		errs() << "}\n";
		
		return;
	}
	
	
	bool StaticSchedulingSamplerCodeGen::runOnFunction(Function &F) {
	
		errs() << " // Start to generating Static Sampling Code (reference based)\n";
		
		/* reading info from previous passes */
		arrayName = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().arrayName;
		arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
		ArrayNameSet = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet;
		LoopRefTNode* LoopRefTree = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().PTLoopRefTree;
		sampleNum = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().sampleNum;
		
		SamplingRate = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().samplingRate;
		
		/* init */
		// initArrayName();
		numberRefToSameArray(LoopRefTree);
		numberLoops(LoopRefTree);
		initIndvName(LoopRefTree);
		
#ifdef PARALLEL
		outloops = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoops;
#endif
		/* generate headers */
		headerGen();

		/* generate rtHistoCal function */
		rtHistoGen();

		/* generate subBlkRT function */
		subBlkRTGen();

#if defined(REFERENCE_GROUP)
		GaussianDistrGen();
		UniformDistrGen();
		rtMergeGen();
#endif
		
#ifdef DumpRTMR
		/* generate rtToMR function */
		rtToMRGen();
		
		/* generate rtDump function */
		rtDumpGen();
		
		/* generate mrDump function */
		mrDumpGen();
#endif
		
		/* generate addr cal function*/
		addrCalFuncGenTop(LoopRefTree);
		
		/* generate rtGen */
		refRTGen(LoopRefTree);
		
		/* generate main function */
		mainGen();

		return false;
	}
	
	void StaticSchedulingSamplerCodeGen::getAnalysisUsage(AnalysisUsage &AU) const {
		AU.setPreservesAll();
		AU.addRequired<idxAnalysis::IndexAnalysis>();
		AU.addRequired<argAnalysis::ArgumentAnalysis>();
		AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
		AU.addRequired<LoopInfoWrapperPass>();
		AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
		AU.addRequired<ivdepAnalysis::IVDependenceAnalysis>();
		AU.addRequired<loopTreeTransform::ParallelLoopTreeTransform>();
		AU.addRequired<AccGraphAnalysis::AccessGraphAnalysis>();
		AU.addRequired<sampleNumAnalysis::SampleNumberAnalysis>();
		return;
	}
	
}
