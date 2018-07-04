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
#include "ssCodeGen_ref.hpp"

#include "llvm/Support/CommandLine.h"

//#define REF_PAIR

using namespace llvm;

static cl::opt<double>
samplingRateEachLoop("sampling rate", cl::Hidden,
                cl::desc("The sampling rate for each loop"));

namespace symFP {
    struct symFP : public FunctionPass {
        static char ID;
        symFP() : FunctionPass(ID) {}
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override {
            
            if (samplingRateEachLoop.getNumOccurrences() > 0) {
                errs() << "sampling rate " << samplingRateEachLoop << "\n";
            }
            
            errs() << " /* Start to analyze function:  \n";
            errs().write_escaped(F.getName()) << " */ \n";
            
            getAnalysis<idxAnalysis::IndexAnalysis>();
            getAnalysis<argAnalysis::ArgumentAnalysis>();
            getAnalysis<gVarAnalysis::GlobalVariableAnalysis>();
            getAnalysis<loopAnalysis::LoopIndvBoundAnalysis>();
            //getAnalysis<brchAnalysis::BranchAnalysis>();

#ifdef REF_PAIR
            getAnalysis<ssCodeGen::StaticSamplingCodeGen>();
#else
            getAnalysis<ssCodeGen_ref::StaticSamplingCodeGen_ref>();
#endif
            
            return false;
        }

        void getAnalysisUsage(AnalysisUsage &AU) const override {
            AU.setPreservesAll();
            AU.addRequired<idxAnalysis::IndexAnalysis>();
            AU.addRequired<argAnalysis::ArgumentAnalysis>();
            AU.addRequired<gVarAnalysis::GlobalVariableAnalysis>();
            AU.addRequired<loopAnalysis::LoopIndvBoundAnalysis>();
            //AU.addRequired<brchAnalysis::BranchAnalysis>();
#ifdef REF_PAIR
            AU.addRequired<ssCodeGen::StaticSamplingCodeGen>();
#else
            AU.addRequired<ssCodeGen_ref::StaticSamplingCodeGen_ref>();
#endif
        }
    };
    
    char symFP::ID = 0;
    static RegisterPass<symFP> X("symFP", "Symbolic Footprint Pass", false, false);
}



