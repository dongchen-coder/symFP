//
//  brchAnalysis.cpp
//  LLVM
//
//  Created by Fangzhou Liu on 2017/7/11.
//
//

#include "brchAnalysis.hpp"

using namespace std;
using namespace llvm;

namespace brchAnalysis {
    char BranchAnalysis::ID = 0;
    static RegisterPass<BranchAnalysis> X("brchAnalysis", "branch condition analysis Pass", false, false);
    
    BranchAnalysis::BranchAnalysis() : FunctionPass(ID) {}
    
    struct LoopBranch {
        Loop *L;
        vector<Value *> BR;
    };
    
    vector<LoopBranch> LoopBranchInfo;
    
    bool BranchAnalysis::runOnFunction(Function &F) {
        
        errs() << "\nStart analysis Branch Conditions \n";
        LoopInfo &LI = getAnalysis<LoopInfoWrapperPass>().getLoopInfo();
        if (!LI.empty()) {
            for(LoopInfo::iterator it = LI.begin(), eit = LI.end(); it != eit; ++it){
                errs() << "Loop:\n";
                
            }
        }

        return false;
    }

    
    void BranchAnalysis::getAnalysisUsage(AnalysisUsage &AU) const {
        AU.setPreservesAll();
        AU.addRequired<LoopInfoWrapperPass>();
        AU.addRequired<idxAnalysis::IndexAnalysis>();
        AU.addRequired<argAnalysis::ArgumentAnalysis>();
        AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
        AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
        return;
    }
}
