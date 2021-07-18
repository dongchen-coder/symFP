//
//  plumCodeGen_ref.cpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 7/13/21.
//

#include "plumCodeGen_ref.hpp"

namespace plumCodeGen_ref {
	char PLUMSamplerCodeGen::ID = 0;
	static RegisterPass<PLUMSamplerCodeGen> X("plumCodeGen_ref", "parallel sampler code generator pass (both reference and iteration). It will choose proper code based on the LIV dependence analysis result", false, false);

	PLUMSamplerCodeGen::PLUMSamplerCodeGen() : FunctionPass(ID) {}

	void PLUMSamplerCodeGen::initIndvName(LoopRefTNode *LoopRefTree) {
		
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

	void PLUMSamplerCodeGen::numberRefToSameArray(LoopRefTNode *LoopRefTree) {
		
		if (LoopRefTree->AA != NULL) {
			if (refToSameArrayCnt.find(LA->arrayName[LoopRefTree->AA]) == refToSameArrayCnt.end()) {
				refNumber[LoopRefTree->AA] = 0;
				refToSameArrayCnt[LA->arrayName[LoopRefTree->AA]] = 1;
			} else {
				refNumber[LoopRefTree->AA] = refToSameArrayCnt[LA->arrayName[LoopRefTree->AA]];
				refToSameArrayCnt[LA->arrayName[LoopRefTree->AA]] += 1;
			}
		}
		
		if (LoopRefTree->next != NULL) {
			for (std::vector<LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
				numberRefToSameArray(*it);
			}
		}
		
		return;
	}

	void PLUMSamplerCodeGen::headerGen() {
		
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
		errs() << "        this->avail_chunk = (trip % chunk_size) == 0 ? trip / chunk_size : trip / chunk_size + 1;\n";
		errs() << "        this->lb = 0;\n";
		errs() << "        this->ub = (this->lb + chunk_size - 1) < this->trip ? (this->lb + chunk_size - 1) : this->trip - 1;\n";
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
		errs() << "        this->ub = (this->lb + chunk_size - 1) < this->trip ? (this->lb + chunk_size - 1) : this->trip - 1;\n";
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
			if (find(visit.begin(), visit.end(), LA->arrayName[it->first]) != visit.end()) { continue; }
			visit.push_back(LA->arrayName[it->first]);
			errs() << "std::unordered_map<uint64_t, uint64_t> LAT_" + LA->arrayName[it->first] + ";\n";
		}
		
		errs() << "std::map<uint64_t, double> RT;\n";
		errs() << "std::map<uint64_t, double> MR;\n";
		return;
	}

	void PLUMSamplerCodeGen::addrCalFuncGenTop(LoopRefTNode *LoopRefTree) {
		std::vector<string> indvs;
		addrCalFuncGen(LoopRefTree, indvs);
		return;
	}

	void PLUMSamplerCodeGen::addrCalFuncGen(LoopRefTNode* LoopRefTree, std::vector<string> indvs) {
		
		if (LoopRefTree->L != NULL) {
			for (std::vector<Value*>::iterator it = LoopRefTree->LIS->IDV->begin(), eit = LoopRefTree->LIS->IDV->end(); it != eit; ++it) {
				indvs.push_back(indvName[(*it)]);
			}
		}
		
		if (LoopRefTree->AA != NULL) {
			errs() << "/* " + LA->arrayName[LoopRefTree->AA] + " " + arrayExpression[LoopRefTree->AA] + " " + std::to_string(refNumber[LoopRefTree->AA])  + " */\n";
			errs() << "int calAddr" + LA->arrayName[LoopRefTree->AA] + std::to_string(refNumber[LoopRefTree->AA]) + "( ";
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

#if 0
	/* Generate the loop iteration vector incremental func */
	void PLUMSamplerCodeGen::LoopIterIncGen(AccGraph *G,
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
	
	void PLUMSamplerCodeGen::refRTGen(LoopRefTNode *LoopRefTree) {
#ifdef PER_ARRAY_ANALYSIS
		std::unordered_set<std::string>::iterator arrayIter = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet.begin();
		for (; arrayIter != getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet.end(); ++arrayIter) {
			errs() << "void interleaving_" << *arrayIter << "() {\n";
			perArrayRTBodyGen(*arrayIter, LoopRefTree);
			errs() << "}\n";
		}
#else
		errs() << "void interleaving() {\n";
		refRTBodyGen(LoopRefTree);
		errs() << "}\n";
#endif
		return;
	}

	void PLUMSamplerCodeGen::mainGen() {
		errs() << "int main() {\n";
		string space = "    ";
		// add timer
		errs() << "#ifdef PAPI_TIMER\n";
		errs() << space << "PAPI_timer_init();\n";
		errs() << space << "PAPI_timer_start();\n";
		errs() << "#endif\n";

#if defined(PER_ARRAY_ANALYSIS)
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
#else
		errs() << space << "interleaving();\n";
#endif
		errs() << "    RTtoMR_AET();\n";
		errs() << "#ifdef PAPI_TIMER\n";
		errs() << "    PAPI_timer_end();\n";
		errs() << "    PAPI_timer_print();\n";
		errs() << "#endif\n";
		errs() << "    rtDump();\n";
		errs() << "    dumpMR();\n";
		
		errs() << "    return 0;\n";
		errs() << "}\n";
		return;
	}

	void PLUMSamplerCodeGen::refRTBodyGen(LoopRefTNode *LoopRefTree) {
		string space = "    ";
		errs() << space << "uint64_t cnt = 0;\n";
		errs() << space << "int tid_to_run, chunk_size;\n";
//		errs() << space << "set<int> candidate_thread_pool;\n";
		errs() << space << "vector<int> candidate_thread_pool, threads_to_exec;\n";
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
			CodeGen->refRTSearchGen(LoopRefTree, LoopRefTree, outloops, "    ");
		} else {
			/* We remove the dummy entry node and run the RTSearchCode for each outermost loop */
			vector<LoopRefTNode*>::iterator loopiter = LoopRefTree->next->begin();
			for (; loopiter != LoopRefTree->next->end(); ++loopiter) {
				LoopRefTNode * subLoop = *loopiter;
				if (subLoop->L && IVA->isParent(subLoop->LIS->IDV->front())) {
					CodeGen->setSchedulingType(DYNAMIC);
					errs() << space << "/* USEING DYNAMIC SCHEUDLING */\n";
				} else {
					CodeGen->setSchedulingType(STATIC);
					errs() << space << "/* USEING STATIC SCHEUDLING */\n";
				}
				CodeGen->refRTSearchGen(subLoop, subLoop, outloops, "    ");
				if (loopiter != (LoopRefTree->next->end() - 1)) {
					/* When change to another outermost loop, we should remove all metadata */
					errs() << space << "candidate_thread_pool.clear();\n";
					errs() << space << "threads_to_exec.clear();\n";
					errs() << space << "progress.clear();\n";
				}
			}
		}
		return;
	}

	void PLUMSamplerCodeGen::perArrayRTBodyGen(std::string array, LoopRefTNode *LoopRefTree) {
		string space = "    ";
		errs() << space << "uint64_t cnt = 0;\n";
		errs() << space << "int tid_to_run, chunk_size;\n";
		errs() << space << "vector<int> candidate_thread_pool, threads_to_exec;\n";
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
			CodeGen->refRTSearchGen(LoopRefTree, LoopRefTree, outloops, "    ");
		} else {
			// We remove the dummy entry node and run the RTSearchCode for each outermost loop
			vector<LoopRefTNode *>::iterator loopiter = LoopRefTree->next->begin();
			for (; loopiter != LoopRefTree->next->end(); ++loopiter) {
				LoopRefTNode * subLoop = *loopiter;
				if (subLoop->L && IVA->isParent(subLoop->LIS->IDV->front())) {
					CodeGen->setSchedulingType(DYNAMIC);
					errs() << space << "/* USEING DYNAMIC SCHEUDLING */\n";
				} else {
					CodeGen->setSchedulingType(STATIC);
					errs() << space << "/* USEING STATIC SCHEUDLING */\n";
				}
				// filter out those outermost loop that does not have the given array access
				if (!LA->DoesLoopContainsArray(subLoop, array))
					continue;
				CodeGen->perArrayRTSearchGen(array, subLoop, subLoop, outloops, "    ");
				if (loopiter != LoopRefTree->next->end()-1) {
					// When change to another outermost loop, we should remove all metadata
					errs() << space << "candidate_thread_pool.clear();\n";
					errs() << space << "threads_to_exec.clear();\n";
					errs() << space << "progress.clear();\n";
				}
			}
		}
		return;
	}
	
	
	bool PLUMSamplerCodeGen::runOnFunction(Function &F) {
	
		errs() << " // Start to generating Static Sampling Code (reference based)\n";
		
		/* reading info from previous passes */
		LA = &getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>();
		IVA = &getAnalysis<ivdepAnalysis::IVDependenceAnalysis>();
		GA = &getAnalysis<AccGraphAnalysis::AccessGraphAnalysis>();
		arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
		ArrayNameSet = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().ArrayNameSet;
		LoopRefTNode* LoopRefTree = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().PTLoopRefTree;
		
		SamplingRate = getAnalysis<sampleNumAnalysis::SampleNumberAnalysis>().samplingRate;
		
		/* init */
		initIndvName(LoopRefTree);
		numberRefToSameArray(LoopRefTree);
		
		CodeGen = new SamplerCodeGenerator(LA, IVA, GA, refNumber);
		if (CodeGen) {
			CodeGen->setSamplingRate(SamplingRate);
		} else {
			return false;
		}
		
#ifdef PARALLEL
		outloops = getAnalysis<loopTreeTransform::ParallelLoopTreeTransform>().outMostLoops;
#endif
		/* generate headers */
		headerGen();

		/* generate rtHistoCal function */
		CodeGen->rtHistoGen();
		/* generate rtToMR function */
		CodeGen->rtToMRGen();
		
		/* generate rtDump function */
		CodeGen->rtDumpGen();
		
		/* generate mrDump function */
		CodeGen->mrDumpGen();
		
		/* generate addr cal function*/
		addrCalFuncGenTop(LoopRefTree);
		
		/* generate rtGen */
		refRTGen(LoopRefTree);
		
		/* generate main function */
		mainGen();

		return false;
	}
	
	void PLUMSamplerCodeGen::getAnalysisUsage(AnalysisUsage &AU) const {
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
