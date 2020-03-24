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
    /* local helper function declaration */
    loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* initThreadNode(int parentLoopLevel);


    loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* initThreadNode(int parentLoopLevel) {
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = (loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *) malloc(sizeof(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode));

        thread_node->isThreadNode = true;
        thread_node->L = NULL;
        thread_node->AA = NULL;
        thread_node->LIS = NULL;
        thread_node->LoopLevel = parentLoopLevel;   
        return thread_node;
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

#if defined(ITER_LEVEL_INTERLEAVING)
    /* 
     * Connect all consecutive reference node with a thread node
     * All thread nodes will be assgined as the child of the last loop node
    */
    void ParallelLoopTreeTransform::LoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator nextIter = LoopRefTree->next->begin();
        // contains the thread node and the loop node
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> * thread_node_vec = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>;
        /* Check if the transform is required
         * if didTransform is false: no extra thread node is added
         * -> no need to switch the LoopRefTree->next
         * else
         * -> have to replace LoopRefTRee->next with thread_node_vec
         */
        while (nextIter != LoopRefTree->next->end()) {
            if ((*nextIter)->AA != NULL) {
                // this is a reference node
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* thread_node_next = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>;
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator refNodeIter = nextIter + 1;
                /* iterate the node after this reference node
                 * if it's also a reference node, append it under the same thread node
                 * until it meets an loop node or reach the end of the children list
                 */
                thread_node_next->push_back(*nextIter);
                while (refNodeIter != LoopRefTree->next->end()) {
                    if ((*refNodeIter)->AA == NULL) { break; }
                    thread_node_next->push_back(*refNodeIter);
                    refNodeIter ++;
                }

                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = initThreadNode(LoopRefTree->LoopLevel);
                thread_node->next = thread_node_next;


                thread_node_vec->push_back(thread_node);

                nextIter = refNodeIter;
            }
            else if ((*nextIter)->L != NULL) { // loop node
                thread_node_vec->push_back(*nextIter);
                nextIter ++;
            }
        }
        // if loop tree has been transformed
        if (!thread_node_vec->empty()) { 
            LoopRefTree->next = thread_node_vec;
        }
    }
#elif defined(ACC_LEVEL_INTERLEAVING) or defined(ALL_LEVEL_INTERLEAVING)
    /* 
     * For each reference node, connect it with a thread node
     * All thread nodes will be assgined as the child of the last loop node
    */
    void ParallelLoopTreeTransform::LoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>::iterator nextIter = LoopRefTree->next->begin();
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* thread_node_vec = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>;
        /* Check if the transform is required
         * if didTransform is false: no extra thread node is added
         * -> no need to switch the LoopRefTree->next
         * else
         * -> have to replace LoopRefTRee->next with thread_node_vec
         */
        bool didTransform = false;
        while (nextIter != LoopRefTree->next->end()) {
            if ((*nextIter)->AA != NULL) { // this is a RefNode
                didTransform = true;
                // create the thread node
                loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = initThreadNode(LoopRefTree->LoopLevel);
                // each thread node connects to a array access node
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *> * thread_node_next = new vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>;
                thread_node_next->push_back(*nextIter);
                thread_node->next = thread_node_next;
                thread_node_vec->push_back(thread_node);
            } else if (didTransform && (*nextIter)->L != NULL) { // this is a LoopNode
                thread_node_vec->push_back(*nextIter);
            }
            nextIter++;
        }
        if (didTransform) {
            LoopRefTree->next = thread_node_vec;
        }
    }
#endif

    void ParallelLoopTreeTransform::tranverseLoopRefTree(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node) {
        if (node != NULL) {
            // this is an loop node, do the transformation
            if (node->L != NULL) {
                errs() << "\t" << node->L->getName() << "\n";
                vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *>* subnodes = node->next;
                LoopTreeTransform(node);
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
