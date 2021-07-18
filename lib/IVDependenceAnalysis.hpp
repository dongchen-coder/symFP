//
//  IVDependenceAnalysis.hpp
//  LLVM
//
//  Created by Fangzhou Liu on 05/06/21.
//
//
#ifndef ivdepAnalysis_hpp
#define ivdepAnalysis_hpp

#include <stdio.h>
#include <string>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/IR/Constants.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "plumCodeGenUtil.hpp"

using namespace llvm;
using namespace std;

/* Debugging flag, enable if debug is needed */
// #define LOOP_DEBUG

namespace ivdepAnalysis {
    struct IVDependenceAnalysis : public FunctionPass {
        static char ID;
        IVDependenceAnalysis();

        SmallSet<DepNode *, 10> DepForest;

        void buildIVDependenceForest(LoopRefTNode * L, SmallSet<DepNode *, 10> & DepForest);
		
		DepNode * getParentOfValue(Value * v, int * type);
		SmallSet<DepNode *, 5> getDependenceOfValue(Value * v);
		
		bool isParent(Value * v);
        bool isLBDep(DepNode * parent, DepNode * child);
        bool isUBDep(DepNode * parent, DepNode * child);
        bool hasDep(DepNode * parent, DepNode * child);
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;

        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}
#endif /* ivdepAnalysis_hpp */
