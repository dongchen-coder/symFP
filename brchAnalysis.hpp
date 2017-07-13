//
//  brchAnalysis.hpp
//  LLVM
//
//  Created by Fangzhou Liu on 2017/7/11.
//
//

#ifndef brchAnalysis_hpp
#define brchAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/CFG.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"

using namespace llvm;

namespace brchAnalysis {
    
    struct BranchAnalysis : public FunctionPass {
        static char ID;
        BranchAnalysis();
        
        /* Extract the Conditions From the Block Within the Loop */
        void BranchConditionAnalysis(Function &F);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        /* Find Branch Conditions in SubLoop */
        void SubLoopAnalysis(Loop *L);
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}


#endif /* brchAnalysis_hpp */
