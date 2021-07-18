//
//  plumSamplerSamplerCodeGenerator.cpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 7/13/21.
//

#include "plumSamplerCodeGenerator.hpp"

void SamplerCodeGenerator::rtHistoGen()
{
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
	return;
}

void SamplerCodeGenerator::rtDumpGen()
{
	
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

void SamplerCodeGenerator::rtToMRGen()
{
	
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

void SamplerCodeGenerator::mrDumpGen()
{
	
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


DepParentInfo SamplerCodeGenerator::getParent(Value * V)
{
	int depType = 0; // 1 means LB dep, -1 means UB dep
	DepNode * parent = IVA->getParentOfValue(V, &depType);
	if (!parent)
		errs() << "/* " << *V << " has no depend parent */\n";
	if (depType < 0)
		errs() << "/* " << *V << " UB depends on " << *(parent->IV->IDV->front()) << " */\n";
	else if (depType > 0)
		errs() << "/* " << *V << " LB depends on " << *(parent->IV->IDV->front()) << " */\n";
	return std::make_pair(parent, depType);
}

/* Generate the loop iteration vector update func */
void SamplerCodeGenerator::LoopIterUpdateGen(AccGraph *G,
									  AccGraphEdge *E,
									  LoopRefTNode* IDomLoop,
									  string space,
									  bool EnableIterationInc,
									  bool EnableIterationUpdate)
{
	// increment part
	if (EnableIterationInc) {
		string upper_bound = this->LA->getBound(IDomLoop->LIS->LB->front().second);
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
			errs() << space << "if (" << progressRepr << " + 1 " << LA->predicateToString(IDomLoop->LIS->PREDICATE->front()) << " " <<  upper_bound << ") {\n";
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
			string lower_bound = LA->getBound((*loopIter)->LIS->LB->front().first);
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

void SamplerCodeGenerator::reuseUpdateGen(LoopRefTNode *access,
								   int currentAccessLoopNestCnt,
								   string accessName,
								   string space)
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
	errs() << space << "    if (LAT_" << LA->arrayName[access->AA] << ".find(access) != LAT_" << LA->arrayName[access->AA] << ".end()) {\n";
	errs() << space << "        rtHistoCal(RT, cnt - LAT_" << LA->arrayName[access->AA] << "[access], 1.0);\n";
	errs() << space << "    }\n";
	errs() << space << "}\n";
#ifdef PARALLEL_CXX_THREAD
	errs() << space << "lck.unlock();\n";
#endif
	errs() << space << "LAT_" << LA->arrayName[access->AA] << "[access] = cnt;\n";
}

void SamplerCodeGenerator::refRTSearchGen(LoopRefTNode *root,
					LoopRefTNode *LoopRefTree,
					std::vector<LoopRefTNode*> & outloops,
					std::string space)
{
	if (this->type == STATIC)
		StaticRefRTSearchGen(root, LoopRefTree, outloops, space);
	else if (this->type == DYNAMIC)
		DynamicRefRTSearchGen(root, LoopRefTree, outloops, space);
	else
		return;
}

void SamplerCodeGenerator::perArrayRTSearchGen(std::string array,
						 LoopRefTNode *root,
						 LoopRefTNode *LoopRefTree,
						 std::vector<LoopRefTNode*> & outloops,
						 std::string space)
{
	if (this->type == STATIC)
		StaticPerArrayRTSearchGen(array, root, LoopRefTree, outloops, space);
	else if (this->type == DYNAMIC)
		DynamicPerArrayRTSearchGen(array, root, LoopRefTree, outloops, space);
	else
		return;
	
}

/* ---------------------
 StatSchSamplerCodeGenerator
 */

/* Generating Reuse Search */
void SamplerCodeGenerator::StaticRefRTSearchGen(LoopRefTNode *root,
												LoopRefTNode *LoopRefTree,
												std::vector<LoopRefTNode*> & outloops,
												std::string space) {
	/**
	 Generate the preparation code at the very beginning of the outermost loop
	 */
	if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
		AccGraph * G = GA->AccGraphLoopMap[LoopRefTree];
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
				// trip counts the number of iterations, default we assume the loop predicate is <
				string trip = LA->getBound(LoopRefTree->LIS->LB->front().second) + " - " + LA->getBound(LoopRefTree->LIS->LB->front().first);
				if (LoopRefTree->LIS->PREDICATE->front() == llvm::CmpInst::ICMP_ULE)
					trip += " + 1";
				errs() << space << "engine = ChunkEngine(chunk_size, (" << trip << "));\n";
				errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
				errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
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
				errs() << space << space << space << space << "Progress p(\"" + LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", {";
				vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
				// traverse the loop from outer to inner
				while (pathRIter != pathToNode.rend()) {
					if (pathRIter == pathToNode.rbegin()) { //
						errs() << "c.first, ";
					} else {
						string lower_bound = LA->getBound((*pathRIter)->LIS->LB->front().first);
						errs() << lower_bound;
						if ( pathRIter != pathToNode.rend()-1)
							errs() << ", ";
					}
					pathRIter++;
				}
				errs() << "}, c);\n";
				errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
				errs() << space << space << space << space << "candidate_thread_pool.erase(candidate_thread_pool.begin());\n";
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
//				errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
				errs() << space << space << space << space << "continue;\n";
				errs() << space << space << space << "}\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
				errs() << "#endif\n";
			}
			string currAccessName = LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
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
			nextAccessName = LA->arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
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
					nextAccessName = LA->arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
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
		errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
		errs() << space << space << space << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
		errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
		errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
		errs() << space << space << space << space << "}\n";
		errs() << space << space << space << "}\n";
		errs() << space << space << "} /* end of thread interleaving loop */\n";
//			errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
		errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
		errs() << space << space << space << "break;\n";
		errs() << space << space << "} /* end of break condition check */\n";
		errs() << space << "} /* end of while(true) */\n";
	} else if (LoopRefTree->AA) { // no loop, just array accesses
		errs() << space << "/* " << LA->arrayName[LoopRefTree->AA] << "[" << LA->arrayExpression[LoopRefTree->AA] << "] */";
//		errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//		errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
//		errs() << space << "}\n";
//		errs() << space << "while (true) {\n";
//#ifdef UNIFORM_INTERLEAVING
//		errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//#endif
//		reuseUpdateGen(LoopRefTree, 0, LA->arrayName[LoopRefTree->AA], space + space);
//		errs() << space << space << "} /* end of thread interleaivng loop */\n";
//		errs() << space << "} /* end of while true */\n";
	}
	return;
}

void SamplerCodeGenerator::StaticPerArrayRTSearchGen(std::string array, LoopRefTNode *root, LoopRefTNode *LoopRefTree, std::vector<LoopRefTNode *> & outloops, std::string space)
{
	/**
	 Generate the preparation code at the very beginning of the outermost loop
	 */
	if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
		AccGraph * G = GA->AccGraphLoopMap[LoopRefTree];
		// If there is no array reference in current LoopNest, we do not need to generate
		// any code.
		vector<LoopRefTNode *> NodesInGraph = G->getAllGraphNodes();
		vector<LoopRefTNode *>::iterator GraphIter = NodesInGraph.begin();
//		vector<LoopRefTNode *>::iterator TargetArrayIter = NodesInGraph.begin();
		LoopRefTNode * currentNode = NULL, * nextNode = NULL;
//		LoopRefTNode * currentTargetArrayNode = NULL, * nextTargetArrayNode = NULL;
		bool reachTheLastAccessNode = false;
		while (true) {
			currentNode = (*GraphIter);
			vector<LoopRefTNode *> pathToNode = G->getPathToNode(currentNode, NULL);
			// generate the preparion work at the very beginning
			// this will be generate if and only if currently we are in the root
			// of the access graph
			if (GraphIter == NodesInGraph.begin()) {
				errs() << space << "engine = ChunkEngine(chunk_size, (" << LA->getBound(LoopRefTree->LIS->LB->front().second) << " - " << LA->getBound(LoopRefTree->LIS->LB->front().first) << "));\n";
				errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
				errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
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
				errs() << space << space << space << space << "Progress p(\"" + LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", {";
		
				vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
				// traverse the loop from outer to inner
				while (pathRIter != pathToNode.rend()) {
					if (pathRIter == pathToNode.rbegin()) { //
						errs() << "c.first, ";
					} else {
						string lower_bound = LA->getBound((*pathRIter)->LIS->LB->front().first);
						errs() << lower_bound;
						if ( pathRIter != pathToNode.rend()-1)
							errs() << ", ";
					}
					pathRIter++;
				}
				errs() << "}, c);\n";
				errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
				errs() << space << space << space << space << "candidate_thread_pool.erase(candidate_thread_pool.begin());\n";
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
//				errs() << space << space << space << space << "candidate_thread_pool.insert(tid_to_run);\n";
				errs() << space << space << space << space << "continue;\n";
				errs() << space << space << space << "}\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
				errs() << "#endif\n";
			}
			string currAccessName = LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
			errs() << space << space << space << "if (progress[tid_to_run].ref == \"" << currAccessName << "\") {\n";
			if (LA->arrayName[currentNode->AA] == array) {
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
			nextAccessName = LA->arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
			
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
					nextAccessName = LA->arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
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
		errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
		errs() << space << space << space << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
		errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
		errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
		errs() << space << space << space << space << "}\n";
		errs() << space << space << space << "}\n";
		errs() << space << space << "} /* end of thread interleaving loop */\n";
//			errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
		errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
		errs() << space << space << space << "break;\n";
		errs() << space << space << "} /* end of break condition check */\n";
		errs() << space << "} /* end of while(true) */\n";
	} else if (LoopRefTree->AA) { // no loop, just array accesses
		errs() << space << "/* " << LA->arrayName[LoopRefTree->AA] << "[" << LA->arrayExpression[LoopRefTree->AA] << "] */";
//		errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//		errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
//		errs() << space << "}\n";
//		errs() << space << "while (true) {\n";
//#ifdef UNIFORM_INTERLEAVING
//		errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//#endif
//		reuseUpdateGen(LoopRefTree, 0, LA->arrayName[LoopRefTree->AA], space + space);
//		errs() << space << space << "} /* end of thread interleaivng loop */\n";
//		errs() << space << "} /* end of while true */\n";
	}
	return;
}





/* ---------------------
 DynScheSamplerCodeGenerator
 */
/* Generating Reuse Search */
void SamplerCodeGenerator::DynamicRefRTSearchGen(
	LoopRefTNode *root,
	LoopRefTNode *LoopRefTree,
	vector<LoopRefTNode*> & outloops,
	string space
	) {
	/**
	 Generate the preparation code at the very beginning of the outermost loop
	 */
	if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
		AccGraph * G = GA->AccGraphLoopMap[LoopRefTree];
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
				errs() << space << "engine = ChunkEngine(chunk_size, (" << LA->getBound(LoopRefTree->LIS->LB->front().second) << " - " << LA->getBound(LoopRefTree->LIS->LB->front().first) << "));\n";
				errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
				errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
				errs() << space << "}\n";
				errs() << space << "while(true) {\n";
#if defined(RANDOM_INTERLEAVING)
				errs() << "chunk_assignment:\n";
#endif
				errs() << space << space << "if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {\n";
				errs() << space << space << space << "int tidx;\n";
				errs() << space << space << space << "while(!candidate_thread_pool.empty()) {\n";
				errs() << space << space << space << space << "tidx = 0;\n";
				errs() << space << space << space << space << "tid_to_run = candidate_thread_pool[tidx];\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] Assign chunk \" << engine.getCurrentChunkRange() << endl;\n";
				errs() << "#endif\n";
				errs() << space << space << space << space << "Chunk c = engine.getNextChunk(tid_to_run);\n";
				errs() << space << space << space << space << "// Init an all-zero vector that equals the size of the iteration vector\n";
				errs() << space << space << space << space << "vector<int> iteration_vector(" << to_string(pathToNode.size()) << ", 0);\n";
				vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
				int iteration_vector_idx = 0;
				// traverse the loop from outer to inner
				while (pathRIter != pathToNode.rend()) {
					if (pathRIter == pathToNode.rbegin()) { //
						errs() << space << space << space << space << "iteration_vector[0] = c.first;\n";
					} else {
						DepParentInfo pinfo = getParent((*pathRIter)->LIS->IDV->front());
						string lower_bound = LA->getBound((*pathRIter)->LIS->LB->front().first);
						if (pinfo.first && pinfo.second > 0) { // has LB dependence
							lower_bound = "iteration_vector[" + to_string(pinfo.first->level - 1) + "]";
						}
						errs() << space << space << space << space << "iteration_vector[" << to_string(iteration_vector_idx) << "] = " << lower_bound << ";\n";
					}
					pathRIter++;
					iteration_vector_idx++;
				}
				errs() << space << space << space << space << "Progress p(\"" + LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", iteration_vector, c);\n";
				errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
				errs() << space << space << space << space << "candidate_thread_pool.erase(candidate_thread_pool.begin() + tidx);\n";
				errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
				errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#endif
				errs() << space << space << space << "} /* end of progress assignment */\n";
				errs() << space << space << "} /* end of chunk availability check */\n";
				/* Generate the Interleaving code for the current outermost loop */
#ifdef UNIFORM_INTERLEAVING
				errs() << space << space << "/* UNIFORM THREAD INTERLEAVING */\n";
				errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//				errs() << space << space << "/* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */\n";
//				errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; });\n";
//				errs() << space << space << "for (auto tid_to_run : threads_to_exec) {\n";
#elif defined(RANDOM_INTERLEAVING)
				errs() << space << space << "/* RANDOM THREAD INTERLEAVING */\n";
				errs() << space << space << "while ( !threads_to_exec.empty()) {\n";
				errs() << space << space << space << "tid_to_run = threads_to_exec[rand() % threads_to_exec.size()];\n";
#endif
				errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] \" << progress[tid_to_run].iteration[0] << \" > \" << progress[tid_to_run].chunk.second << endl;\n";
				errs() << "#endif\n";
				errs() << space << space << space << space << "continue;\n";
				errs() << space << space << space << "}\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
				errs() << "#endif\n";
			}
			string currAccessName = LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
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
			nextAccessName = LA->arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
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
				for (auto edge : edges) { // iterate the edges from inner to outer
					carryLoop = edge->getCarryLoop();
					sink = edge->getSinkNode();
					nextAccessName = LA->arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
					if (carryLoop) { // CASE 1 or CASE 3
						// CASE 3
						if (G->QueryEngine.getImmediateLoopDominator(currentNode) == carryLoop) {
							// increment the LIV only
							errs() << "/* CASE 3 */\n";
							LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, !G->QueryEngine.areTwoAccessInSameLevel(currentNode, sink));
							errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
							errs() << space << space << space << space << space << "continue;\n";
							errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
						} else { // CASE 1
							errs() << "/* CASE 1 */\n";
							// pop iteration vector L1.size() times
							// push LB of loop in L2 in iteration vector
							// and increment the LIV
							LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, true);
							errs() << space << space << space << space << space << "progress[tid_to_run].increment(\"" << nextAccessName << "\");\n";
							errs() << space << space << space << space << space << "continue;\n";
							errs() << space << space << space << space << "} /* end of check to " << currAccessName << " */\n";
						}
					} else if (edge->hasSource(currentNode) && edge->hasSink(nextNode)) { // CASE 2
						errs() << "/* CASE 2 */\n";
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
		errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
		errs() << "#ifdef SIMULATOR_DEBUG\n";
		errs() << space << space << space << space << space << "cout << \"Put \" << tid_to_run << \" in candidate_thread_pool\" << endl;\n";
		errs() << "#endif\n";
		errs() << space << space << space << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
		errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
		errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
		errs() << space << space << space << space << "}\n";
		errs() << space << space << space << "}\n";
		errs() << space << space << "} /* end of thread interleaving loop */\n";
//		errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
		errs() << space << space << "cout << candidate_thread_pool.size() << endl;\n";
		errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
		errs() << space << space << space << "break;\n";
		errs() << space << space << "} /* end of break condition check */\n";
		errs() << space << "} /* end of while(true) */\n";
	} else if (LoopRefTree->AA) {
		reuseUpdateGen(LoopRefTree, 0, LA->arrayName[LoopRefTree->AA], space);
	}
	return;
}

void SamplerCodeGenerator::DynamicPerArrayRTSearchGen(string array,
	LoopRefTNode *root,
	LoopRefTNode *LoopRefTree,
	vector<LoopRefTNode*> & outloops,
	string space
	) {
	/**
	 Generate the preparation code at the very beginning of the outermost loop
	 */
	if (find(outloops.begin(), outloops.end(), LoopRefTree) != outloops.end()) {
		AccGraph * G = GA->AccGraphLoopMap[LoopRefTree];
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
				errs() << space << "engine = ChunkEngine(chunk_size, (" << LA->getBound(LoopRefTree->LIS->LB->front().second) << " - " << LA->getBound(LoopRefTree->LIS->LB->front().first) << "));\n";
				errs() << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
				errs() << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
				errs() << space << "}\n";
				errs() << space << "while(true) {\n";
#if defined(RANDOM_INTERLEAVING)
				errs() << "chunk_assignment:\n";
#endif
				errs() << space << space << "if (!candidate_thread_pool.empty() && engine.hasNextChunk()) {\n";
				errs() << space << space << space << "int tidx;\n";
				errs() << space << space << space << "while(!candidate_thread_pool.empty()) {\n";
				errs() << space << space << space << space << "tidx = 0;\n";
				errs() << space << space << space << space << "tid_to_run = candidate_thread_pool[tidx];\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] Assign chunk \" << engine.getCurrentChunkRange() << endl;\n";
				errs() << "#endif\n";
				errs() << space << space << space << space << "Chunk c = engine.getNextChunk(tid_to_run);\n";
				errs() << space << space << space << space << "// Init an all-zero vector that equals the size of the iteration vector\n";
				errs() << space << space << space << space << "vector<int> iteration_vector(" << to_string(pathToNode.size()) << ", 0);\n";
				vector<LoopRefTNode *>::reverse_iterator pathRIter = pathToNode.rbegin();
				int iteration_vector_idx = 0;
				// traverse the loop from outer to inner
				while (pathRIter != pathToNode.rend()) {
					if (pathRIter == pathToNode.rbegin()) { //
						errs() << space << space << space << space << "iteration_vector[0] = c.first;\n";
					} else {
						DepParentInfo pinfo = getParent((*pathRIter)->LIS->IDV->front());
						string lower_bound = LA->getBound((*pathRIter)->LIS->LB->front().first);
						if (pinfo.first && pinfo.second > 0) { // has LB dependence
							lower_bound = "iteration_vector[" + to_string(pinfo.first->level - 1) + "]";
						}
						errs() << space << space << space << space << "iteration_vector[" << to_string(iteration_vector_idx) << "] = " << lower_bound << ";\n";
					}
					pathRIter++;
					iteration_vector_idx++;
				}
				errs() << space << space << space << space << "Progress p(\"" + LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]) << "\", iteration_vector, c);\n";
				errs() << space << space << space << space << "progress[tid_to_run] = p;\n";
				errs() << space << space << space << space << "candidate_thread_pool.erase(candidate_thread_pool.begin() + tidx);\n";
				errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
				errs() << space << space << space << space << "threads_to_exec.emplace_back(tid_to_run);\n";
#endif
				errs() << space << space << space << "} /* end of progress assignment */\n";
				errs() << space << space << "} /* end of chunk availability check */\n";
				/* Generate the Interleaving code for the current outermost loop */
#ifdef UNIFORM_INTERLEAVING
				errs() << space << space << "/* UNIFORM THREAD INTERLEAVING */\n";
				errs() << space << space << "for (tid_to_run = 0; tid_to_run < THREAD_NUM; tid_to_run++) {\n";
//				errs() << space << space << "/* RANDOMLY CHOOSE THREAD TO RUN EACH REFERENCE */\n";
//				errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; });\n";
//				errs() << space << space << "for (auto tid_to_run : threads_to_exec) {\n";
#elif defined(RANDOM_INTERLEAVING)
				errs() << space << space << "/* RANDOM THREAD INTERLEAVING */\n";
				errs() << space << space << "while ( !threads_to_exec.empty()) {\n";
				errs() << space << space << space << "tid_to_run = threads_to_exec[rand() % threads_to_exec.size()];\n";
#endif
				errs() << space << space << space << "if (!progress[tid_to_run].isInBound()) {\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << space << "cout << \"[\" << tid_to_run << \"] \" << progress[tid_to_run].iteration[0] << \" > \" << progress[tid_to_run].chunk.second << endl;\n";
				errs() << "#endif\n";
				errs() << space << space << space << space << "continue;\n";
				errs() << space << space << space << "}\n";
				errs() << "#ifdef SIMULATOR_DEBUG\n";
				errs() << space << space << space << "cout << \"[\" << tid_to_run << \"] Iterate \" << progress[tid_to_run].ref << \" at \" << progress[tid_to_run].getIteration() << endl;\n";
				errs() << "#endif\n";
			}
			string currAccessName = LA->arrayName[currentNode->AA] + to_string(refNumber[currentNode->AA]), nextAccessName = "";
			errs() << space << space << space << "if (progress[tid_to_run].ref == \"" << currAccessName << "\") {\n";
			if (LA->arrayName[currentNode->AA] == array) {
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
			nextAccessName = LA->arrayName[nextNode->AA] + to_string(refNumber[nextNode->AA]);
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
					nextAccessName = LA->arrayName[edge->getSinkNode()->AA] + to_string(refNumber[edge->getSinkNode()->AA]);
					if (carryLoop) { // CASE 1 or CASE 3
						// CASE 3
						if (G->QueryEngine.getImmediateLoopDominator(currentNode) == carryLoop) {
							// increment the LIV only
							LoopIterUpdateGen(G, edge, carryLoop, space + space + space + space, true, !G->QueryEngine.areTwoAccessInSameLevel(currentNode, sink));
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
		errs() << space << space << space << space << "if (find(candidate_thread_pool.begin(), candidate_thread_pool.end(), tid_to_run) == candidate_thread_pool.end()) {\n";
		errs() << space << space << space << space << space << "candidate_thread_pool.emplace_back(tid_to_run);\n";
#ifdef RANDOM_INTERLEAVING
		errs() << space << space << space << space << space << "threads_to_exec.erase(find(threads_to_exec.begin(), threads_to_exec.end(), tid_to_run));\n";
		errs() << space << space << space << space << space << "goto chunk_assignment;\n";
#endif
		errs() << space << space << space << space << "}\n";
		errs() << space << space << space << "}\n";
		errs() << space << space << "} /* end of thread interleaving loop */\n";
//		errs() << space << space << "random_shuffle(threads_to_exec.begin(), threads_to_exec.end(), [](int n) { return rand() % n; } );\n";
		errs() << space << space << "cout << candidate_thread_pool.size() << endl;\n";
		errs() << space << space << "if (candidate_thread_pool.size() == THREAD_NUM && !engine.hasNextChunk()) {\n";
		errs() << space << space << space << "break;\n";
		errs() << space << space << "} /* end of break condition check */\n";
		errs() << space << "} /* end of while(true) */\n";
	}
	return;
}
