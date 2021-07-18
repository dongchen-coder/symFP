//
//  IVDependenceAnalysis.cpp
//  LLVM
//
//  Created by Fangzhou Liu on 05/06/21.
//
//
#include <stack>
#include <unordered_set>
#include "llvm/Support/raw_ostream.h"
#include "IVDependenceAnalysis.hpp"


#define DEBUG_TYPE "ivdependence"

using namespace llvm;
using namespace std;


typedef unordered_set<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> USet;
typedef stack<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> Stack;
namespace ivdepAnalysis {


	IVDependenceAnalysis::IVDependenceAnalysis() : FunctionPass(ID) {}

	DepNode * buildIVDependenceNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * L)
	{
		return new DepNode(L->LoopLevel, L->LIS);
	}

	void IVDependenceAnalysis::buildIVDependenceForest(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * L, SmallSet<DepNode *, 10> & DepForest)
	{
		/* remove the dummy loop node */
		while (!L->LIS) {
			L = L->next->front();
		}
		USet visited;
		Stack LoopStack;
		stack<DepNode *> NodeStack;
		DepNode * Root = NULL, * tmpHead = NULL;
		DenseMap<DepNode *, DepNode *> ParentMap;
		if (L && !L->AA) {
			LoopStack.push(L);
			Root = buildIVDependenceNode(L);
			NodeStack.push(Root);
		}
		while (!LoopStack.empty()) {
			L = LoopStack.top();
			tmpHead = NodeStack.top();
			LoopStack.pop();
			NodeStack.pop();
			if (!L->AA && visited.find(L) == visited.end()) {
				visited.insert(L);
			}
			for (auto nextL : *(L->next)) {
				if (!nextL->AA && visited.find(nextL) == visited.end()) {
					LoopStack.push(nextL);
					DepNode * nextNode = buildIVDependenceNode(nextL);
					DepNode * tmpChild = nextNode;
					NodeStack.push(nextNode);
					ParentMap.insert(make_pair(nextNode, tmpHead));
					while (ParentMap.find(tmpChild) != ParentMap.end()) {
						DepNode * parent = ParentMap[tmpChild];
						if (hasDep(parent, nextNode)/* L and nextL has dependence */) {
							parent->insertDependence(nextNode);
							nextNode->insertParent(parent);
							if (DepForest.find(parent) == DepForest.end())
								DepForest.insert(parent);
						}
						tmpChild = parent;
					}
				}
			}
		}
	}

	bool IVDependenceAnalysis::isParent(Value * v)
	{
		for (auto node : DepForest) {
			if (node->IV->IDV->front() == v) {
				return !node->Dependences.empty();
			}
		}
		return false;
	}

	bool IVDependenceAnalysis::isLBDep(DepNode * parent, DepNode * child)
	{
		if (parent->IV)
//		// k = i+1
			LLVM_DEBUG(dbgs() << "LB Parent: " << *(parent->IV->IDV->front()) << "\n");
		if (child->IV)
			LLVM_DEBUG(dbgs() << "LB Child: " << *(child->IV->LB->front().first) << "\n");
//		return false;
		return (getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(child->IV->LB->front().first).find(getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(parent->IV->IDV->front())) != string::npos);

	}

	bool IVDependenceAnalysis::isUBDep(DepNode * parent, DepNode * child)
	{
		// k <= i
		if (parent->IV)
			LLVM_DEBUG(dbgs() << "UB Parent: " << *(parent->IV->IDV->front()) << "\n");
		if (child->IV)
			LLVM_DEBUG(dbgs() << "UB Child: " << *(child->IV->LB->front().second) << "\n");
		return (getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(child->IV->LB->front().second).find(getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(parent->IV->IDV->front())) != string::npos);
	}

	bool IVDependenceAnalysis::hasDep(DepNode * a, DepNode * b)
	{
		LLVM_DEBUG(
				   dbgs() << "check dependence " << *(a->IV->IDV->front())
				   << " with " << *(b->IV->IDV->front()) << "\n";
		);
		if (a->level == b->level)
			return false;
		DepNode * parent = NULL, * child = NULL;
		if (a->level > b->level) {
			parent = b;
			child = a;
		} else {
			parent = a;
			child = b;
		}
		return isLBDep(parent, child) || isUBDep(parent, child);
	}

	SmallSet<DepNode *, 5> IVDependenceAnalysis::getDependenceOfValue(Value * v)
	{
		for (auto node : DepForest) {
			if (node->IV->IDV->front() == v) {
				return node->Dependences;
			}
		}
		return SmallSet<DepNode *, 5>();
	}

	DepNode * IVDependenceAnalysis::getParentOfValue(Value * v, int * type)
	{
		for (auto node : DepForest) {
			for (auto dep : node->Dependences) {
				if (dep->IV->IDV->front() == v) {
					if (isLBDep(node, dep)) {
						* type = 1;
					} else if (isUBDep(node, dep)) {
						* type = -1;
					} else {
						*type = 0;
					}
					return node;
				}
			}
		}
		return NULL;
		
	}

	/* Main function for loop induction variable analysis and loop bound analysis
	 * This pass is to construct a tree structure that stores loop hierarchy and references
	 */
	bool IVDependenceAnalysis::runOnFunction(Function &F) 
	{
		LoopRefTNode * Root = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
		if (!Root->L && !Root->AA) {
			for (auto subLoop : *(Root->next)) {
				if (subLoop->L)
					buildIVDependenceForest(subLoop, DepForest);
			}
		}
		// dump all dependence
		errs() << "/* Start IV Dependence Analysis ";
		for (auto node : DepForest) {
			errs() << getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(node->IV->IDV->front()) << ":\n";
			for (auto subNode : node->Dependences) {
				errs() << "  - "
				<< getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().getBound(subNode->IV->IDV->front())
				<< "(" << subNode->level << ")\n";
			}
		}
		errs() << "Finish to analyze IV Dependence */";
		return false;
	}

	 /* Dependence relation of the analysis paths */
	void IVDependenceAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
		AU.setPreservesAll();
		AU.addRequired<idxAnalysis::IndexAnalysis>();
		AU.addRequired<argAnalysis::ArgumentAnalysis>();
		AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
		AU.addRequired<LoopInfoWrapperPass>();
		AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
		return;
	}
	
	char IVDependenceAnalysis::ID = 0;
	static RegisterPass<IVDependenceAnalysis> X("ivdepAnalysis", "induction variable dependence analysis Pass", false, false);


}
