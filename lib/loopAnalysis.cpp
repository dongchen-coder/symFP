//
//  loopAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include "loopAnalysis.hpp"

using namespace llvm;
using namespace std;

namespace loopAnalysis {
    
    char LoopIndvBoundAnalysis::ID = 0;
    static RegisterPass<LoopIndvBoundAnalysis> X("loopAnalysis", "loop induction variable/bound analysis Pass", false, false);
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    
    string LoopIndvBoundAnalysis::predicateToString(llvm::CmpInst::Predicate PREDICATE) {
        switch (PREDICATE) {
            case llvm::CmpInst::ICMP_EQ:
                return "=";
                break;
            case llvm::CmpInst::ICMP_NE:
                return "!=";
                break;
            case llvm::CmpInst::ICMP_SGE:
                return ">=";
                break;
            case llvm::CmpInst::ICMP_SGT:
                return ">";
                break;
            case llvm::CmpInst::ICMP_SLE:
                return "<=";
                break;
            case llvm::CmpInst::ICMP_SLT:
                return "<";
                break;
            case llvm::CmpInst::ICMP_UGE:
                return ">=";
                break;
            case llvm::CmpInst::ICMP_UGT:
                return ">";
                break;
            case llvm::CmpInst::ICMP_ULE:
                return "<=";
            case llvm::CmpInst::ICMP_ULT:
                return ">";
            default:
                break;
        }
        return "";
    }
    
    LoopIndvBoundAnalysis::LoopBound LoopIndvBoundAnalysis::findLoopBound(Loop *L, Value *var) {
        
        Value* ub;
        Value* lb;
        
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        
        for (BasicBlock *b : L->getBlocks()) {
            if ((b->getName().str().compare(0, 8, "for.cond") == 0) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                
                for(BasicBlock::iterator it = b->begin(), eit = b->end(); it != eit; ++it) {
                    if (isa<CmpInst>(*it)) {
                        
                        Value *v = it->getOperand(1);
                        if (isa<ConstantInt>(v)) {
                            ub = v;
                        } else if (isa<LoadInst>(v)) {
                            LoadInst* ldTmp = dyn_cast<LoadInst>(v);
                            ub = ldTmp->getOperand(0);
                        } else {
                            // the upper bound is an expression
                            ub = v;
                        }
                    }
                }
                
                // iterate all its current predecessors
                for(auto pit = pred_begin(b), pet = pred_end(b); pit != pet; ++pit) {
                    BasicBlock *pred = *pit;
                    if (!(pred->getName().str().compare(0, 7, "for.inc") == 0)) {
                        // iterate all store instructions, find the first operand when the second operand is a idv
                        for (BasicBlock::iterator it = pred->begin(), eit = pred->end(); it != eit; ++it) {
                            if (isa<StoreInst>(*it)) {
                                Value *v = it->getOperand(0);
                                if (it->getOperand(1)->getName().equals(var->getName())) {
//                                    errs() << "Lower Bound for Loop is: " << dyn_cast<ConstantInt>(v)->getValue() <<"\n";
                                    lb = v;
                                }
                            }
                        }
                    }
                }
            }
            
        }
        
        return make_pair(lb, ub);
    }

    /* Get expression for bounds */
    string LoopIndvBoundAnalysis::getBound(Value *bound) {
        
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
                    break;
                default:
                    break;
            }
            
        } else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        
        return "";
    }
    
    /* Get all condition basic block among all basic blocks inside loop L */
    vector<BasicBlock *> LoopIndvBoundAnalysis::getSubLoopCondBlock(Loop *L) {
        vector<BasicBlock *> temp;
        for (Loop *SL : L->getSubLoops()) {
            for (BasicBlock *BB: SL->getBlocks()) {
                //if (std::regex_match (BB->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                if (BB->getName().str().compare(0, 8, "for.cond") == 0) {
                    temp.push_back(BB);
                }
            }
        }
        return temp;
    }
    
    /* Get all increment basic block among all basic blocks inside loop L */
    vector<BasicBlock *> LoopIndvBoundAnalysis::getSubLoopIncBlock(Loop *L) {
        vector<BasicBlock *> temp;
        for (Loop *SL : L->getSubLoops()) {
            for (BasicBlock *BB: SL->getBlocks()) {
                if (BB->getName().str().compare(0, 7, "for.inc") == 0) {
                    temp.push_back(BB);
                }
            }
        }
        return temp;
    }
    
    /* Create and init node for references */
    LoopIndvBoundAnalysis::LoopRefTNode* LoopIndvBoundAnalysis::ExtractRefInfo(Instruction* Inst) {
        
        LoopIndvBoundAnalysis::LoopRefTNode* node = NULL;
        
        if (isa<LoadInst>(Inst))  {
            LoadInst* LD = dyn_cast<LoadInst>(Inst);
            if (isa<GetElementPtrInst>(LD->getOperand(0))) {
                node = (LoopRefTNode*) malloc(sizeof(LoopRefTNode));
                node->isThreadNode = false;
                node->L = NULL;
                node->AA = &(*Inst);
                node->LIS = NULL;
                node->LoopLevel = node->LoopLevel+1;
                node->next = NULL;
            }
        }
        if (isa<StoreInst>(Inst)) {
            StoreInst* ST = dyn_cast<StoreInst>(Inst);
            if (isa<GetElementPtrInst>(ST->getOperand(1))) {
                node = (LoopRefTNode*) malloc(sizeof(LoopRefTNode));
                node->isThreadNode = false;
                node->L = NULL;
                node->AA = &(*Inst);
                node->LIS = NULL;
                node->LoopLevel = node->LoopLevel+1;
                node->next = NULL;
            }
        }
        
        return node;
    }
    
    /* Decroate the loop tree with references */
    LoopIndvBoundAnalysis::LoopRefTNode* LoopIndvBoundAnalysis::LoopTreeConstructionRef(LoopRefTNode* root, vector<BasicBlock*> BBList) {
#ifdef LOOP_DEBUG
        errs() << "\n start to analysis Ref\n";
#endif
        vector<LoopRefTNode*>* nextWithRef = new vector<LoopRefTNode*>;
        
        if (root->next == NULL) {

#ifdef LOOP_DEBUG
            errs() << "Last level loop\n";
#endif
            for (vector<BasicBlock*>::iterator BB = BBList.begin(), eBB = BBList.end(); BB != eBB; ++BB) {
                for (BasicBlock::iterator Inst = (*BB)->begin(), eInst = (*BB)->end(); Inst != eInst; ++Inst) {
                    LoopRefTNode* Tmp = ExtractRefInfo(&*Inst);
                    if (Tmp != NULL) {
                        nextWithRef->push_back(Tmp);
                    }
                }
            }
            
        } else {
            
#ifdef LOOP_DEBUG
            errs() << "not Last level loop\n";
#endif
            
            vector<LoopRefTNode*>::iterator nextIter = root->next->begin();
            vector<BasicBlock*> BBinLoop = (*nextIter)->L->getBlocks();

#ifdef LOOP_DEBUG
            if (root->L != NULL) {
                root->L->dump();
            }
            errs() << "--";
            (*nextIter)->L->dump();
#endif

            for (vector<BasicBlock*>::iterator it = BBList.begin(), eit = BBList.end(); it != eit; ++it) {
            
                BasicBlock* BB = *it;

#ifdef LOOP_DEBUG
                errs() << BB->getName() << "\n";
#endif
                
                if (find(BBinLoop.begin(), BBinLoop.end(), BB) == BBinLoop.end() || BBinLoop.empty()) {
                    
                    for (BasicBlock::iterator Inst = BB->begin(), eInst = BB->end(); Inst != eInst; ++Inst) {

                        LoopRefTNode* Tmp = ExtractRefInfo(&*Inst);
                        if (Tmp != NULL) {
                            nextWithRef->push_back(Tmp);
                        }
                    }
                    
                } else {
                    
                    BBinLoop.erase(find(BBinLoop.begin(), BBinLoop.end(), BB));
                    if (BBinLoop.empty()) {
                        
                        nextWithRef->push_back(*nextIter);
                        ++nextIter;
                        if (nextIter != root->next->end()) {
                            BBinLoop = (*nextIter)->L->getBlocks();
                        }
                    }
                    
                }
            }
        }
        
        if (root->next != NULL) {
            delete(root->next);
        }
        root->next = nextWithRef;
        
        /* add Ref nodes to sub loops */
        for (vector<LoopRefTNode*>::iterator it = root->next->begin(), eit = root->next->end(); it != eit; ++it) {
            if ((*it)->L != NULL) {
                
                vector<BasicBlock*> newBBList;
                for (vector<BasicBlock*>::iterator bit = BBList.begin(), ebit = BBList.end(); bit != ebit; ++bit) {
                    if (find((*it)->L->getBlocks(), *bit) != (*it)->L->getBlocks().end()) {
                        newBBList.push_back(*bit);
                    }
                }
                
                (*it) = LoopTreeConstructionRef((*it), newBBList);
            }
        }
        
        return root;
    }
    
    /* Init loopInfoStruct for each loop L */
    LoopIndvBoundAnalysis::LoopInfoStruct* LoopIndvBoundAnalysis::ExtractLoopInfo(Loop *L) {
        
        vector<Value *>* IDV = new vector<Value *>;
        vector<LoopIndvBoundAnalysis::LoopBound>* LB = new vector<LoopIndvBoundAnalysis::LoopBound>;
        vector<Value *>* INC = new vector<Value *>;
        vector<llvm::CmpInst::Predicate>* PREDICATE = new vector<llvm::CmpInst::Predicate>;
        
        LoopIndvBoundAnalysis::LoopInfoStruct* LIS = new LoopIndvBoundAnalysis::LoopInfoStruct;
        
        /* Extracting loop bounds (LB) / induction variable (IDV) / stride (INC) */
        vector<BasicBlock *> subLoopIncBlocks = getSubLoopIncBlock(L);
        for (BasicBlock *b : L->getBlocks()) {
            if ((b->getName().str().compare(0, 7, "for.inc") == 0) && find(subLoopIncBlocks.begin(), subLoopIncBlocks.end(), b) == subLoopIncBlocks.end()) {
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        IDV->push_back(II.getOperand(0));
                        LB->push_back(findLoopBound(L, II.getOperand(0)));
                    }
                    
                    if (isa<StoreInst>(II)) {
                        INC->push_back(II.getOperand(0));
                    }
                }
            }
        }
        
        /* Extracting loop condition predicate */
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        for (BasicBlock *b : L->getBlocks()) {
            if ((b->getName().str().compare(0, 8, "for.cond") == 0) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                
                PREDICATE->push_back(llvm::CmpInst::BAD_ICMP_PREDICATE);
                
                for(BasicBlock::iterator it = b->begin(), eit = b->end(); it != eit; ++it) {
                    if (isa<CmpInst>(*it)) {
                        CmpInst* CmpTmp = dyn_cast<CmpInst>(it);
                        PREDICATE->pop_back();
                        PREDICATE->push_back(CmpTmp->getPredicate());
                    }
                }
            }
        }
        
        /* assign extracted information to struct */
        LIS->IDV = IDV;
        LIS->LB = LB;
        LIS->INC = INC;
        LIS->PREDICATE = PREDICATE;
        
        return LIS;
    }
    
    /* Construct loop tree sub routine (references are not included): extract information for sub loops */
    void LoopIndvBoundAnalysis::LoopTreeConstruction(Loop* L, LoopRefTNode * root, int level) {
        
        if (L != NULL) {
            
            root->LoopLevel = level;
            root->isThreadNode = false;
            root->L = L;
            root->LIS = ExtractLoopInfo(L);
            root->AA = NULL;
            if (L->getSubLoopsVector().empty()) {
                root->next = NULL;
            } else {
                root->next = new vector<LoopRefTNode *>;
            }
            
            for (Loop *SL : L->getSubLoops()) {
                LoopRefTNode * subTmp = (LoopRefTNode *) malloc(sizeof(LoopRefTNode));
                root->next->push_back(subTmp);
                LoopTreeConstruction(SL, subTmp, level+1);
            }
        }
        
        return;
    }
    
    
    /* Construct loop tree (references are not included) */
    LoopIndvBoundAnalysis::LoopRefTNode* LoopIndvBoundAnalysis::LoopTreeConstructionLoop(LoopRefTNode* root) {
        
        /* using llvm internal loop analysis */
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        
        if (!LI.empty()) {
            root = (LoopRefTNode *) malloc(sizeof(LoopRefTNode));
            root->LoopLevel = 0;
            root->isThreadNode = false;
            root->L = NULL;
            root->LIS = NULL;
            root->AA = NULL;
            root->next = new vector<LoopRefTNode *>;
            
            for(LoopInfo::reverse_iterator it = LI.rbegin(), eit = LI.rend(); it != eit; ++it){
                
                LoopRefTNode * Tmp = (LoopRefTNode *) malloc(sizeof(LoopRefTNode));
                Tmp->L = (*it);
                Tmp->isThreadNode = false;

#ifdef LOOP_DEBUG
                errs() << "\n";
                Tmp->L->dump();
                errs() << "\n";
#endif
                
                Tmp->LIS = ExtractLoopInfo(*it);
                Tmp->AA = NULL;
                Tmp->LoopLevel = 1;
                if (Tmp->L->getSubLoops().empty()) {
                    Tmp->next = NULL;
                } else {
                    Tmp->next = new vector<LoopRefTNode *>;
                }
                
                root->next->push_back(Tmp);
                
                for (Loop *SL : (*it)->getSubLoops()) {
                    LoopRefTNode * subTmp = (LoopRefTNode *) malloc(sizeof(LoopRefTNode));
                    Tmp->next->push_back(subTmp);
                    LoopTreeConstruction(SL, subTmp, 2);
                }
            }
        }
        
        return root;
    }
    
    /* Dump loop/reference tree */
    void LoopIndvBoundAnalysis::DumpLoopTree(LoopRefTNode* LTroot, std::string prefix) {
   
#ifdef LOOP_DEBUG
        errs() << "Start to dump loop tree\n";
#endif
        
        if (LTroot != NULL) {
            
            if (LTroot->LIS != NULL) {
                for (vector<Value*>::iterator it = LTroot->LIS->IDV->begin(), eit = LTroot->LIS->IDV->end(); it != eit; ++it ) {
                    if ((*it) != NULL) {
                        errs() << prefix << (*it)->getName() << "\n";
                    } else {
                        errs() << prefix << "NULL idv\n";
                    }
                }
                
                for (vector<LoopBound>::iterator it = LTroot->LIS->LB->begin(), eit = LTroot->LIS->LB->end(); it != eit; ++it) {
                    errs() << prefix << "Loop Bound: (" << getBound((*it).first) << ", " << getBound((*it).second) << ")\n";
                }
                
                for (vector<Value*>::iterator it = LTroot->LIS->INC->begin(), eit = LTroot->LIS->INC->end(); it != eit; ++it) {
                    errs() << prefix << "Loop inc: " << getBound(*it) << "\n";
                }
                
                for (vector<llvm::CmpInst::Predicate>::iterator it = LTroot->LIS->PREDICATE->begin(), eit = LTroot->LIS->PREDICATE->end(); it != eit; ++it) {
                    errs() << prefix << "Loop predicate: " << predicateToString(*it) << "\n";
                }
            }
            
            if (LTroot->AA != NULL) {
                errs() << prefix + "array access ";
                errs() << arrayName[LTroot->AA] << " ";
                errs() << arrayExpression[LTroot->AA] << "\n";
            }
            
            if (LTroot->next != NULL) {
                for (vector<LoopRefTNode*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
                    DumpLoopTree(*it, prefix + "--");
                }
            }
        } else {
            errs() << "LTroot == NULL\n";
        }
        
        return;
    }
    
    /* Return all the basic blocks in a function */
    vector<BasicBlock*> getBasicBlocks(Function &F) {
        vector<BasicBlock*> BBList;
        
        for (BasicBlock &BB : F) {
            BBList.push_back(&BB);
        }
        
        return BBList;
    }
    
    /* Main function for loop induction variable analysis and loop bound analysis
     * This pass is to construct a tree structure that stores loop hierarchy and references
     */
    bool LoopIndvBoundAnalysis::runOnFunction(Function &F) {
        
        errs() << "\n /* Start analysis loops\n";
        
        arrayName = getAnalysis<idxAnalysis::IndexAnalysis>().arrayName;
        arrayExpression = getAnalysis<idxAnalysis::IndexAnalysis>().arrayExpression;
        
        LoopRefTNode* LTroot = NULL;
        
        /* construct loop tree */
        LTroot = LoopTreeConstructionLoop(LTroot);
        
        /* decroate the loop tree with references */
        LTroot = LoopTreeConstructionRef(LTroot, getBasicBlocks(F));

        LoopRefTree = LTroot;
        
        /* dump loopRef tree */
        DumpLoopTree(LTroot, "");
        
        errs() << "\nFinish analysis loops */ \n";
        
        return false;
    }
    
    /* Dependence relation of the analysis paths */
    void LoopIndvBoundAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<LoopInfoWrapperPass>();
        return;
    }
    
}
