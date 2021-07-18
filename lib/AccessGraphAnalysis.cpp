//
//  AccessGraphAnalysis.cpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 5/10/21.
//
#include <stack>
#include <unordered_set>
#include "llvm/Support/raw_ostream.h"
#include "AccessGraphAnalysis.hpp"

#define DEBUG_TYPE "access-graph-analysis"

using namespace llvm;
using namespace std;
namespace AccGraphAnalysis {

	AccessGraphAnalysis::AccessGraphAnalysis() : FunctionPass(ID) {}
	void AccessGraphAnalysis::buildAccessGraph(LoopRefTNode * L)
	{
		if (!L)
			return;
		if (!L->L && !L->AA) {
			for (auto outermostLoop : *(L->next)) {
				if (outermostLoop->L) {
					AccGraphLoopMap[outermostLoop] = new AccGraph(outermostLoop);
				}
			}
		}
	}

	void AccessGraphAnalysis::dumpAccessGraph(AccGraph * G, string prefix)
	{
		for (auto edge: G->getAllGraphEdges()) {
			errs() << getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().arrayName[edge->getSourceNode()->AA]
			<< " -> "
			<<getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().arrayName[edge->getSinkNode()->AA]
			<< "\n";
		}
	}


	/* Main function for loop induction variable analysis and loop bound analysis
	 * This pass is to construct a tree structure that stores loop hierarchy and references
	 */
	bool AccessGraphAnalysis::runOnFunction(Function &F)
	{
		buildAccessGraph(getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree);
		unordered_map<LoopRefTNode *, AccGraph *>::iterator graphIter = AccGraphLoopMap.begin();
		errs() << "/* Start to analysis the access graph";
		while (graphIter != AccGraphLoopMap.end()) {
			dumpAccessGraph(graphIter->second, "");
			graphIter ++;
			errs() << "\n";
		}
		errs() << "*/\n";
		return false;
	}

	 /* Dependence relation of the analysis paths */
	void AccessGraphAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
		AU.setPreservesAll();
		AU.addRequired<idxAnalysis::IndexAnalysis>();
		AU.addRequired<argAnalysis::ArgumentAnalysis>();
		AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
		AU.addRequired<LoopInfoWrapperPass>();
		AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
		return;
	}

	char AccessGraphAnalysis::ID = 0;
	static RegisterPass<AccessGraphAnalysis> X("lagAnalysis", "access graph analysis Pass", false, false);
}
