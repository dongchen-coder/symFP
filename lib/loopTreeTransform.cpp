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

    /* compute loop iteration space */
    uint64_t computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);

    string getBound(Value *bound);


    loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* initThreadNode(int parentLoopLevel) {
        
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode * thread_node = (loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *) malloc(sizeof(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode));

        thread_node->isThreadNode = true;
        thread_node->L = NULL;
        thread_node->AA = NULL;
        thread_node->LIS = NULL;
        thread_node->LoopLevel = parentLoopLevel;   
        return thread_node;
    }

    string getLoopInc(Value *inc) {
        if (isa<Instruction>(inc)) {
            Instruction *inst = cast<Instruction>(inc);
            switch (inst->getOpcode()) {
                case Instruction::Add:
                case Instruction::Sub:
                case Instruction::Mul:
                case Instruction::FDiv:
                case Instruction::SDiv:
                case Instruction::UDiv:
                    if (isa<ConstantInt>(inst->getOperand(0))) {
                        return getLoopInc(inst->getOperand(0));
                    } else {
                        return getLoopInc(inst->getOperand(1));
                    }
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(inc)) {
            return to_string(dyn_cast<ConstantInt>(inc)->getValue().getSExtValue());
        }
        return "";
    }

    string getBound(Value *bound) {
        
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
                    return inst->getName();
                    break;
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }

    uint64_t computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree) {
        uint64_t iterSpace = 0;
        if (LoopRefTree->next != NULL) {
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = LoopRefTree->next->begin(), eit = LoopRefTree->next->end(); it != eit; ++it) {
                if ((*it)->L == NULL) {
                    iterSpace += 1;
                } else {
                    iterSpace += computeIterSpace(*it);
                }
            }
            iterSpace *= ((uint64_t) stoi(getBound((*LoopRefTree->LIS->LB)[0].second)) - (uint64_t) stoi(getBound((*LoopRefTree->LIS->LB)[0].first))) / (uint64_t) stoi(getLoopInc((*LoopRefTree->LIS->INC)[0]));        
        }
        return iterSpace;
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

    void iterateArrayInstr(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node, vector<Instruction*> & arrayVec) {
        if (node->L == NULL) {
            arrayVec.push_back(node->AA);
            return;
        } else if (node->next != NULL) {
            for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator lit = node->next->begin(); lit != node->next->end(); ++lit) {
                iterateArrayInstr(*lit, arrayVec);
            }
        }
    }

    void ParallelLoopTreeTransform::computePerIterationSpace() {
        for (vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>::iterator it = outMostLoops.begin(), eit = outMostLoops.end(); it != eit; ++it) {
            loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* omloop = *it;
            // go over this outermost loop
            // compute the iteration space for each reference within this outermost loop
            uint64_t IS = computeIterSpace(*it);
            if (omloop->L != NULL) {
                vector<Instruction*> arrayVec;
                iterateArrayInstr(omloop, arrayVec);
                for(vector<Instruction*>::iterator ait = arrayVec.begin(); ait != arrayVec.end(); ++ait) {
                    outMostLoopPerIterationSpace[*ait] = IS;
                }
            }
        }
        return;
    }

#if defined(ITER_LEVEL_INTERLEAVING)
    /* 
     * Connect all consecutive reference node with a thread node
     * All thread nodes will be assigned as the child of the last loop node
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
        computePerIterationSpace();

#if defined(ACC_LEVEL_INTERLEAVING) || defined(ITER_LEVEL_INTERLEAVING)
        errs() << "\n /* Start transform loop tree\n";
        
        tranverseLoopRefTree(PTLoopRefTree);
// #ifdef DEBUG
        DumpLoopTree(PTLoopRefTree, "");
// #endif
        errs() << "\nFinish transform loop tree */ \n";
#endif

        
        
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
