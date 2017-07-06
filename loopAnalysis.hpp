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
#include "llvm/IR/Constants.h"


#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"

using namespace llvm;
using namespace std;

namespace loopAnalysis {
    struct LoopIndvBoundAnalysis : public FunctionPass {
        static char ID;
        LoopIndvBoundAnalysis();
        
        void subLoop(Loop *L);
        /* Find all Basic Induction Variable */
        void findIDV(Loop *L);
        void findLoopBound(Loop *L);
        void inductionVariableAnalysis(Function &F);
        vector<BasicBlock *> getSubLoopCondBlock(Loop *L);
        void dumpLoopIDV();
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
        
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}



#endif /* loopAnalysis_hpp */
