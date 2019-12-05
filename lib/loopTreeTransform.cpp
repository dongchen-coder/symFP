//
//  loopTreeTransform.cpp
//  LLVM
//
//  Created by Fangzhou Liu on 12/03/19.
//
//

#include "loopTreeTransform.hpp"

using namespace llvm;
using namespace std;

namespace loopTreeTransform {
    
    char ParallelLoopTreeTransform::ID = 0;
    static RegisterPass<ParallelLoopTreeTransform> X("loopTreeTransform", "loop tree transform for parallel program Pass", false, false);
    
    ParallelLoopTreeTransform::ParallelLoopTreeTransform() : FunctionPass(ID) {}

    bool containsAANode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node);
    bool isLastLevelLoopNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node);

    /* Check if the node contains an Array Access Node */
    bool containsAANode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        for(vector<TreeNodeBase *>::iterator it = node->next->begin(); it != node->next->end(); ++it) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* next_node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
            if (next_node != nullptr && next_node->AA != NULL) {
                return true;
            }
        }
        return false;
    }

    /* Check if the node is the last level loop node */
    bool isLastLevelLoopNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        for(vector<TreeNodeBase *>::iterator it = node->next->begin(); it != node->next->end(); ++it) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* next_node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
            if (next_node != nullptr && next_node->L != NULL) {
                return false;
            }
        }
        return true;
    }

    /* Get all the out loops */
    /* The thread local iteration calculation code
        will be generated based on all out loops */
   void ParallelLoopTreeTransform::findAllOutMostLOops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        for (vector<TreeNodeBase*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(*it);
            if (node != nullptr && node->L != NULL) {
                outMostLoops.push_back(node);
            }
        }
#ifdef DEBUG
        errs() << "/* # of Out-most Loops: " << outMostLoops.size() << " */ \n";  
#endif
        return;
    }

    void ParallelLoopTreeTransform::insertThreadNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {

        vector<TreeNodeBase *>* tmp_next = LoopRefTree->next;
        vector<TreeNodeBase *>* thread_node_vec = new vector<TreeNodeBase *>;

        ParallelLoopTreeTransform::ThreadNode* thread_node = (ThreadNode*) malloc(sizeof(ThreadNode));
    
        thread_node->next = tmp_next;
        thread_node_vec->push_back(thread_node);
        LoopRefTree->next = thread_node_vec;

    }

#if defined(RANDOM_INTERLEAVING)
    void ParallelLoopTreeTransform::RandomLoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {

    }
#else
    /* 
     * Call when reach the last level loop node: loop node whose childen are all reference node 
     * For each reference node, connect it with a thread node
     * All thread nodes will be assgined as the child of the last loop node
    */
    void ParallelLoopTreeTransform::LoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        vector<TreeNodeBase *>::iterator nextIter = LoopRefTree->next->begin();
        vector<TreeNodeBase *>* thread_node_vec = new vector<TreeNodeBase*>;
        while (nextIter != LoopRefTree->next->end()) {
            // create the thread node
            ParallelLoopTreeTransform::ThreadNode* thread_node = (ThreadNode*) malloc(sizeof(ThreadNode));
            // each thread node connects to a array access node
            vector<TreeNodeBase *> * thread_node_next = new vector<TreeNodeBase *>;
            thread_node_next->push_back(*nextIter);
            thread_node->next = thread_node_next;
            thread_node_vec->push_back(thread_node);
        }
        LoopRefTree->next = thread_node_vec;
    }
#endif

    void ParallelLoopTreeTransform::tranverseLoopRefTree(TreeNodeBase* node) {
        if (node != NULL) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * lnd = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(node);
            // this is an loop node, do the transformation
            if (lnd != nullptr && lnd->L != NULL) {
                if (lnd->L != NULL) {
                    errs() << lnd->L->getName();
#if defined(UNIFORM_INTERLEAVING)
                    if (containsAANode(lnd)) {
                        insertThreadNode(lnd);
                    }
#elif defined(RANDOM_INTERLEAVING)

#else               
                    if (isLastLevelLoopNode(lnd)) {
                        LoopTreeTransform(lnd); // insert thread node for each reference node 
                    } else if (containsAANode(lnd)) {
                        insertThreadNode(lnd); // insert thread node only for the reference node
                    } 

#endif
                    vector<TreeNodeBase *>* subnodes = lnd->next;
                    for (vector<TreeNodeBase *>::iterator nit=subnodes->begin(), neit=subnodes->end(); nit!=neit; ++nit) {
                        tranverseLoopRefTree(*nit);
                    }
                }
            }
        }
    }


    /* Dump loop/reference tree */
    void ParallelLoopTreeTransform::DumpLoopTree(TreeNodeBase* LTroot, string prefix) {
   
#ifdef LOOP_DEBUG
        errs() << "Start to dump loop tree\n";
#endif
        if (LTroot != NULL) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node = static_cast<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>(LTroot);
            if (node != nullptr) {    
                if (node->LIS != NULL) {
                    errs() << prefix << "--------------" << "\n";
                    errs() << prefix << "|  LoopNode  |" << "\n";
                    errs() << prefix << "--------------" << "\n";
                }
            
                if (node->AA != NULL) {
                    errs() << prefix << "--------------" << "\n";
                    errs() << prefix << "| AccessNode |" << "\n";
                    errs() << prefix << "--------------" << "\n";
                }
            } else {
                ThreadNode* node = static_cast<ThreadNode*>(LTroot);
                if (node != nullptr) {
                    errs() << prefix << "--------------" << "\n";
                    errs() << prefix << "| ThreadNode |" << "\n";
                    errs() << prefix << "--------------" << "\n";
                }
            }
            
            if (LTroot->next != NULL) {
                for (vector<TreeNodeBase*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
                    DumpLoopTree(*it, prefix + "--");
                }
            }
        } else {
            errs() << "LTroot == NULL\n";
        }
        
        return;
    }
    
    
    /* Main function for loop induction variable analysis and loop bound analysis
     * This pass is to construct a tree structure that stores loop hierarchy and references
     */
    bool ParallelLoopTreeTransform::runOnFunction(Function &F) {

        PTLoopRefTree = getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>().LoopRefTree;
        
        errs() << "\n /* Start transform loop tree\n";
        
        tranverseLoopRefTree(PTLoopRefTree);

// #ifdef DEBUG
        DumpLoopTree(PTLoopRefTree, "");
// #endif
        
        errs() << "\nFinish transform loop tree */ \n";
        
        return false;
    }
    
    /* Dependence relation of the analysis paths */
    void ParallelLoopTreeTransform::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        return;
    }
    
}
