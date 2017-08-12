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
    
    
    
    LoopIndvBoundAnalysis::LoopBound LoopIndvBoundAnalysis::findLoopBound(Loop *L, Value *var) {
        Value* ub;
        Value* lb;
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match(b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                for(BasicBlock::iterator it = b->begin(), eit = b->end(); it != eit; ++it) {
                    if (isa<ICmpInst>(*it)) {
//                        errs() << "Upper Bound for Loop is: ";
                        Value *v = it->getOperand(1);
                        if (isa<ConstantInt>(v)) {
//                            errs() << dyn_cast<ConstantInt>(v)->getValue();
                            ub = v;
                        }
                    }
                }
                // iterate all its current predecessors
                for(auto pit = pred_begin(b), pet = pred_end(b); pit != pet; ++pit) {
                    BasicBlock *pred = *pit;
                    if (!std::regex_match(pred->getName().str(), std::regex("^for.inc\\d*$|^for.inc$"))) {
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
                default:
                    break;
            }
        }
        else if (isa<ConstantInt>(bound)) {
            return to_string(dyn_cast<ConstantInt>(bound)->getValue().getSExtValue());
        }
        return "";
    }

    
    
    
    /*
     * 获得Loop L所有子Loop中的cond BasicBlock
     */
    
    vector<BasicBlock *> LoopIndvBoundAnalysis::getSubLoopCondBlock(Loop *L) {
        vector<BasicBlock *> temp;
        for (Loop *SL : L->getSubLoops()) {
            for (BasicBlock *BB: SL->getBlocks()) {
                if (std::regex_match (BB->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                    //                    errs() << BB->getName() + " \n";
                    temp.push_back(BB);
                }
            }
        }
        return temp;
    }
    
    LoopIndvBoundAnalysis::LoopTreeNodes* LoopIndvBoundAnalysis::ExtractRefInfo(Instruction* Inst) {
        
        LoopIndvBoundAnalysis::LoopTreeNodes* node = NULL;
        
        if (isa<LoadInst>(Inst))  {
            LoadInst* LD = dyn_cast<LoadInst>(Inst);
            if (isa<GetElementPtrInst>(LD->getOperand(0))) {
                node = (LoopTreeNodes*) malloc(sizeof(LoopTreeNodes));
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
                node = (LoopTreeNodes*) malloc(sizeof(LoopTreeNodes));
                node->L = NULL;
                node->AA = &(*Inst);
                node->LIS = NULL;
                node->LoopLevel = node->LoopLevel+1;
                node->next = NULL;
            }
        }
        
        return node;
    }
    
    LoopIndvBoundAnalysis::LoopTreeNodes* LoopIndvBoundAnalysis::LoopTreeConstructionRef(LoopTreeNodes* root, vector<BasicBlock*> BBList) {
        
        vector<LoopTreeNodes*>* nextWithRef = new vector<LoopTreeNodes*>;
        
        if (root->next == NULL) {
            
            for (vector<BasicBlock*>::iterator BB = BBList.begin(), eBB = BBList.end(); BB != eBB; ++BB) {
                for (BasicBlock::iterator Inst = (*BB)->begin(), eInst = (*BB)->end(); Inst != eInst; ++Inst) {
                    LoopTreeNodes* Tmp = ExtractRefInfo(&*Inst);
                    if (Tmp != NULL) {
                        nextWithRef->push_back(Tmp);
                    }
                }
            }
            
        } else {
            
            vector<LoopTreeNodes*>::iterator nextIter = root->next->begin();
            vector<BasicBlock*> BBinLoop = (*nextIter)->L->getBlocks();
            
            for (vector<BasicBlock*>::iterator it = BBList.begin(), eit = BBList.end(); it != eit; ++it) {
            
                BasicBlock* BB = *it;
                
                if (find(BBinLoop.begin(), BBinLoop.end(), BB) == BBinLoop.end() || BBinLoop.empty()) {
                    
                    for (BasicBlock::iterator Inst = BB->begin(), eInst = BB->end(); Inst != eInst; ++Inst) {

                        LoopTreeNodes* Tmp = ExtractRefInfo(&*Inst);
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
        for (vector<LoopTreeNodes*>::iterator it = root->next->begin(), eit = root->next->end(); it != eit; ++it) {
            if ((*it)->L != NULL) {
                (*it) = LoopTreeConstructionRef((*it), (*it)->L->getBlocks());
            }
        }
        
        return root;
    }
    
    
    LoopIndvBoundAnalysis::LoopInfoStruct* LoopIndvBoundAnalysis::ExtractLoopInfo(Loop *L) {
        
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        
        vector<Value *>* IDV = new vector<Value *>;
        vector<LoopIndvBoundAnalysis::LoopBound>* LB = new vector<LoopIndvBoundAnalysis::LoopBound>;
        
        LoopIndvBoundAnalysis::LoopInfoStruct* LIS = new LoopIndvBoundAnalysis::LoopInfoStruct;
        
        // find all BasicBlocks in the loop L
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match (b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)")) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        IDV->push_back(II.getOperand(0));
                        LB->push_back(findLoopBound(L, II.getOperand(0)));
                    }
                }
            }
        }
        
        LIS->IDV = IDV;
        LIS->LB = LB;
        
        return LIS;
    }
    
    void LoopIndvBoundAnalysis::LoopTreeConstruction(Loop* L, LoopTreeNodes * root, int level) {
        
        if (L != NULL) {
            
            //root = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
            root->LoopLevel = level;
            root->L = L;
            root->LIS = ExtractLoopInfo(L);
            root->AA = NULL;
            if (L->getSubLoopsVector().empty()) {
                root->next = NULL;
            } else {
                root->next = new vector<LoopTreeNodes *>;
            }
            
            for (Loop *SL : L->getSubLoops()) {
                LoopTreeNodes * subTmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                root->next->push_back(subTmp);
                LoopTreeConstruction(SL, subTmp, level+1);
            }
        }
        
        return;
    }
    
    LoopIndvBoundAnalysis::LoopTreeNodes* LoopIndvBoundAnalysis::LoopTreeConstructionLoop(LoopTreeNodes* root) {
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        
        if (!LI.empty()) {
            root = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
            root->LoopLevel = 0;
            root->L = NULL;
            root->LIS = NULL;
            root->AA = NULL;
            root->next = new vector<LoopTreeNodes *>;
            
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                
                LoopTreeNodes * Tmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                Tmp->L = (*it);
                Tmp->LIS = ExtractLoopInfo(*it);
                Tmp->AA = NULL;
                Tmp->LoopLevel = 1;
                Tmp->next = new vector<LoopTreeNodes *>;
                
                root->next->push_back(Tmp);
                
                for (Loop *SL : (*it)->getSubLoops()) {
                    LoopTreeNodes * subTmp = (LoopTreeNodes *) malloc(sizeof(LoopTreeNodes));
                    Tmp->next->push_back(subTmp);
                    LoopTreeConstruction(SL, subTmp, 2);
                }
            }
        }
        
        return root;
    }
    
    /* Dump loopRef tree */
    void LoopIndvBoundAnalysis::DumpLoopTree(LoopTreeNodes* LTroot, std::string prefix) {
        
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
            }
            
            if (LTroot->AA != NULL) {
                errs() << prefix + "array access ";
                LTroot->AA->dump();
            }
            
            if (LTroot->next != NULL) {
                for (vector<LoopTreeNodes*>::iterator it = LTroot->next->begin(), eit = LTroot->next->end(); it != eit; ++it) {
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
     * 
     * This pass is to construct a tree structure that stores loop hierarchy and references
     *
     */
    bool LoopIndvBoundAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis loops\n";
        
        LoopTreeNodes* LTroot = NULL;
        
        /* construct loop tree */
        LTroot = LoopTreeConstructionLoop(LTroot);

        /* decroate the loop tree with references */
        LTroot = LoopTreeConstructionRef(LTroot, getBasicBlocks(F));
        
        DumpLoopTree(LTroot, "");
        
        errs() << "\nFinish analysis loops\n";
        
        return false;
    }
    
    /* Dependence relation of the analysis paths */
    void LoopIndvBoundAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        return;
    }
    
}
