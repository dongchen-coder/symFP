//
//  loopAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/19/17.
//
//

#ifndef loopAnalysis_hpp
#define loopAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/CFG.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"

using namespace llvm;

namespace loopAnalysis {
    struct LoopIndvBoundAnalysis : public FunctionPass {
        static char ID;
        LoopIndvBoundAnalysis();
        
        void subLoop(Loop *L);
        void inductionVariableAnalysis(Function &F);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* loopAnalysis_hpp */
