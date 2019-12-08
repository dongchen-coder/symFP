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
        for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = node->next->begin(); it != node->next->end(); ++it) {
            if ((*it)->AA != NULL) {
                return true;
            }
        }
        return false;
    }

    /* Check if the node is the last level loop node */
    bool isLastLevelLoopNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        for(vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator it = node->next->begin(); it != node->next->end(); ++it) {
            if ((*it)->L != NULL) {
                return false;
            }
        }
        return true;
    }

    /* Get all the out loops */
    /* The thread local iteration calculation code
        will be generated based on all out loops */
   void ParallelLoopTreeTransform::findAllOutMostLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
            if ((*it)->L != NULL) {
                outMostLoops.push_back((*it));
            }
        }
#ifdef DEBUG
        errs() << "/* # of Out-most Loops: " << outMostLoops.size() << " */ \n";  
#endif
        return;
    }

    void ParallelLoopTreeTransform::insertThreadNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {

        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* thread_node_vec = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>;

        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = (loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *) malloc(sizeof(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode));

        thread_node->isThreadNode = true;
        thread_node->L = NULL;
        thread_node->AA = NULL;
        thread_node->LIS = NULL;
        thread_node->LoopLevel = -1;
        thread_node->next = LoopRefTree->next;
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
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator nextIter = LoopRefTree->next->begin();
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* thread_node_vec = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>;
        while (nextIter != LoopRefTree->next->end()) {
            // this is a RefNode
            if ((*nextIter)->AA != NULL) {

            }
            // create the thread node
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = (loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *) malloc(sizeof(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode));

            thread_node->isThreadNode = true;
            thread_node->L = NULL;
            thread_node->AA = NULL;
            thread_node->LIS = NULL;
            thread_node->LoopLevel = -1;
            // each thread node connects to a array access node
            vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> * thread_node_next = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>;
            thread_node_next->push_back(*nextIter);
            thread_node->next = thread_node_next;
            thread_node_vec->push_back(thread_node);
            nextIter++;
        }
        LoopRefTree->next = thread_node_vec;
    }
#endif

    void ParallelLoopTreeTransform::tranverseLoopRefTree(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        if (node != NULL) {
            // this is an loop node, do the transformation
            if (node->L != NULL) {
                errs() << node->L->getName();
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* subnodes = node->next;
#if defined(UNIFORM_INTERLEAVING)
                if (containsAANode(node)) {
                    insertThreadNode(node);
                }
#elif defined(RANDOM_INTERLEAVING)

#else               
                if (isLastLevelLoopNode(node)) {
                    errs() << "isLastLevelLoopNode\n"; 
                    LoopTreeTransform(node); // insert thread node for each reference node 
                } else if (containsAANode(node)) {
                    errs() << "containsAANode\n"; 
                    insertThreadNode(node); // insert thread node only for the reference node
                } 

#endif
                for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator nit=subnodes->begin(), neit=subnodes->end(); nit!=neit; ++nit) {
                    tranverseLoopRefTree(*nit);
                }
            } else if (node->next != NULL) {
                for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator nit=node->next->begin(), neit=node->next->end(); nit!=neit; ++nit) {
                    tranverseLoopRefTree(*nit);
                }
            }
        }
    }


    /* Dump loop/reference tree */
    void ParallelLoopTreeTransform::DumpLoopTree(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LTroot, string prefix) {
   
#ifdef LOOP_DEBUG
        errs() << "Start to dump loop tree\n";
#endif
        if (LTroot != NULL) {
            if (LTroot->isThreadNode) {    
                errs() << prefix << "--------------" << "\n";
                errs() << prefix << "| ThreadNode |" << "\n";
                errs() << prefix << "--------------" << "\n";
            } else {
                if (LTroot->LIS != NULL) {
                    errs() << prefix << "--------------" << "\n";
                    errs() << prefix << "|  LoopNode  |" << "\n";
                    errs() << prefix << "--------------" << "\n";
                }
            
                if (LTroot->AA != NULL) {
                    errs() << prefix << "--------------" << "\n";
                    errs() << prefix << "| AccessNode |" << "\n";
                    errs() << prefix << "--------------" << "\n";
                }
            }
            
            if (LTroot->next != NULL) {
                for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
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
        findAllOutMostLoops(PTLoopRefTree);
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