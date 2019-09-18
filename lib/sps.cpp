#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Instructions.h"
#include "llvm/IR/User.h"
#include "llvm/IR/IntrinsicInst.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/IR/Constants.h"
#include <string>
#include <vector>
#include <set>
#include <regex>

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
//#include "brchAnalysis.hpp"
//#include "ssCodeGen.hpp"

/* sequential static sampling gen */
#include "ssCodeGen_ref.hpp"
/* sequential resverior sampling gen */
#include "rsCodeGen_ref.hpp"
/* parallel 1:1 interleaving static sampling gen */
#include "psCodeGen_ref.hpp"


//#define REF_PAIR

using namespace llvm;

namespace symFP {
    struct symFP : public FunctionPass {
        static char ID;
        symFP() : FunctionPass(ID) {}
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override {
            
            errs() << " /* Analyze function: ";
            errs().write_escaped(F.getName()) << " */ \n";
            
            return false;
        }

        void getAnalysisUsage(AnalysisUsage &AU) const override {
            AU.setPreservesAll();
#if defined(RS_SAMPLE)
            AU.addRequired<rsCodeGen_ref::ReservoirSamplingCodeGen_ref>();
#elif defined(PARALLEL)
            AU.addRequired<psCodeGen_ref::ParallelSamplingCodeGen_ref>();
#elif defined(REF_PAIR)
            AU.addRequired<ssCodeGen::StaticSamplingCodeGen>();
#else
            AU.addRequired<ssCodeGen_ref::StaticSamplingCodeGen_ref>();
#endif
        }
    };
    
    char symFP::ID = 0;
    static RegisterPass<symFP> X("sps", "Static Parallel Sampling Pass", false, false);
}



