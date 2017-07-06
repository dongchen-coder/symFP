//
//  loopAnalysis.cpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#include <regex>
#include "loopAnalysis.hpp"

using namespace llvm;
using namespace std;

namespace loopAnalysis {
    
    char LoopIndvBoundAnalysis::ID = 0;
    static RegisterPass<LoopIndvBoundAnalysis> X("loopAnalysis", "loop induction variable/bound analysis Pass", false, false);
    
    typedef pair<StringRef, vector<Value *>> LoopIDVInfo;
    
    vector<LoopIDVInfo> LoopIDVInfoList;
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    int loopcont = 0;
    void LoopIndvBoundAnalysis::subLoop(Loop *L) {
        for (Loop *SL : L->getSubLoops()) {
            loopcont ++;
            errs() << "Sub Loop:\n";
            findIDV(SL);
            subLoop(SL);
        }
        
        return;
    }
    
    void LoopIndvBoundAnalysis::findIDV(Loop *L) {
        vector<BasicBlock *> subLoopCondBlocks = getSubLoopCondBlock(L);
        errs() << "BasicBlocks in SubLoop is: ";
        for (BasicBlock *BB : subLoopCondBlocks) {
            errs() << BB->getName() + " ";
        }
        vector<Value *> temp;
        // find all BasicBlocks in the loop L
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match (b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)")) && find(subLoopCondBlocks.begin(), subLoopCondBlocks.end(), b) == subLoopCondBlocks.end()) {
                errs() << b->getName() + " \n";
                
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        errs() << "IDV for Loop " + L->getName().str() + " is " + II.getOperand(0)->getName() + "\n";
                        temp.push_back(II.getOperand(0));
                    }
                }
            }
        }
        LoopIDVInfoList.push_back(make_pair(L->getName(), temp));
        errs() << "\n";
    }
    
    void LoopIndvBoundAnalysis::findLoopBound(Loop *L) {
        for (BasicBlock *b : L->getBlocks()) {
            if (std::regex_match(b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                for(BasicBlock::iterator it = b->begin(), eit = b->end(); it != eit; ++it) {
                    if (isa<ICmpInst>(*it)) {
                        errs() << "Upper Bound for Loop is: ";
                        Value *v = it->getOperand(1);
                        if (isa<ConstantInt>(v))
                        {
                            Value *use = it->getOperand(0);
                            errs() << *(dyn_cast<Instruction>(use)->getOperand(0));
                            Instruction *def = dyn_cast<Instruction>(use);
                            
                            errs() << "def: " << *def << "\n";
                            for(Value::use_iterator i = def->use_begin(), ie = def->use_end(); i!=ie; ++i) {
                                Value *v = *i;
                                Instruction *vi = dyn_cast<Instruction>(v);
                                errs() << "\t\t" << *vi << "\n";
                            }
//                            errs() << *(dyn_cast<Instruction>(def)->getUse);
//                            errs() << cast<ConstantInt>(v)->getValue();
                        }
                    }
                }
            }
        }
    }
    
    void LoopIndvBoundAnalysis::inductionVariableAnalysis(Function &F) {
        errs() << "Induction variable analysis\n";
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        if (!LI.empty()) {
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                errs() << "Loop:\n";
                findIDV(*it);
                findLoopBound(*it);
                subLoop(*it);
                loopcont ++;
            }
        }
        dumpLoopIDV();
        errs() << "There's ";
        errs() << loopcont;
        errs() << " Loops In total\n";
        return;
    }
    
    bool LoopIndvBoundAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis loops\n";
        
        inductionVariableAnalysis(F);
        
        return false;
    }
    
    void LoopIndvBoundAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        return;
    }
    
    /*
     * 获得Loop L所有子Loop中的cond BasicBlock
     */
    vector<BasicBlock *> LoopIndvBoundAnalysis::getSubLoopCondBlock(Loop *L) {
        vector<BasicBlock *> temp;
        for (Loop *SL : L->getSubLoops()) {
            for (BasicBlock *BB: SL->getBlocks()) {
                if (std::regex_match (BB->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                    errs() << BB->getName() + " \n";
                    temp.push_back(BB);
                }

            }
        }
        return temp;
    }
    
    void LoopIndvBoundAnalysis::dumpLoopIDV() {
        for (LoopIDVInfo info: LoopIDVInfoList) {
            errs() << "For Loop " + info.first.str() + ": ";
            for (Value *v : info.second) {
                errs() << v->getName() + " ";
            }
            errs() << "\n";
        }
    }
}
