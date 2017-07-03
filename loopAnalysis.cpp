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
    
    
    LoopIndvBoundAnalysis::LoopIndvBoundAnalysis() : FunctionPass(ID) {}
    int loopcont = 0;
    void LoopIndvBoundAnalysis::subLoop(Loop *L) {
        for (Loop *SL : L->getSubLoops()) {
            loopcont ++;
            findIDV(SL);
            errs() << "Loop:\n";
            PHINode* phi = SL->getCanonicalInductionVariable();
            if (phi != NULL) {
                phi->dump();
            } else {
                errs() << "  Can not find induction variable\n";
            }
            subLoop(SL);
        }
        
        return;
    }
    
    void LoopIndvBoundAnalysis::findIDV(Loop *L) {
        
        vector<BasicBlock*> loopbase = L->getBlocks();
        // find all BasicBlocks in the loop L
//        errs() << "Loop " + L->getName() + " base Basic Blocks are ";
        for (BasicBlock * b : loopbase) {
            if (std::regex_match (b->getName().str(), std::regex("^for.cond$|^for.cond\\d*$)"))) {
                errs() << b->getName() + " \n";
                for (Instruction &II : *b) {
                    if (isa<LoadInst>(II)) {
                        errs() << "IDV is " + II.getOperand(0)->getName() + "\n";
                    }
                }
            }
        }
        errs() << "\n";
    }
    
    void LoopIndvBoundAnalysis::inductionVariableAnalysis(Function &F) {
        errs() << "Induction variable analysis\n";
        
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        if (!LI.empty()) {
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                errs() << "Loop:\n";
                findIDV(*it);
                PHINode* phi = (*it)->getCanonicalInductionVariable();
                if (phi != NULL) {
                    phi->dump();
                } else {
                    errs() << " Can not find induction variable\n";
                }
                subLoop(*it);
                loopcont ++;
            }
        }
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
}
