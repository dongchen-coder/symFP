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
// #include "rsCodeGen_ref.hpp"
#if defined(PARALLEL)
#   include "loopTreeTransform.hpp"
#   if defined(UNIFORM_INTERLEAVING)
        /* parallel 1:1 interleaving static sampling gen */
#       include "uiAccCodeGen_ref.hpp"
#       include "uiIterCodeGen_ref.hpp"
#       include "uiAccCodeGenOpt_ref.hpp"
#   elif defined(RANDOM_INTERLEAVING)
        /* parallel random interleaving static sampling gen */
#       include "riAccCodeGen_ref.hpp"
#       include "riIterCodeGen_ref.hpp"
#		include "riCodeGen_ref.hpp"
#   elif defined(THREAD_TUNING)
#       include "staticTuningCodeGen_ref.hpp"
#   else 
#       include "modelCodeGen_ref.hpp"
#       include "modelCodeGenOpt_ref.hpp"
#   endif
#endif


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
#if defined(PARALLEL) 
            AU.addRequired<loopTreeTransform::ParallelLoopTreeTransform>();
#   if defined(UNIFORM_INTERLEAVING) && defined(ITER_LEVEL_INTERLEAVING)
            AU.addRequired<uiIterCodeGen_ref::IterLevelUISamplingCodeGen_ref>();
#   elif defined(UNIFORM_INTERLEAVING) && defined(ACC_LEVEL_INTERLEAVING)
            AU.addRequired<uiAccCodeGenOpt_ref::AccLevelUISamplingCodeGenOpt_ref>();
            // AU.addRequired<uiAccCodeGen_ref::AccLevelUISamplingCodeGen_ref>();
#   elif defined(RANDOM_INTERLEAVING) && defined(ITER_LEVEL_INTERLEAVING)
            AU.addRequired<riIterCodeGen_ref::IterLevelRICodeGen_ref>();
#   elif defined(RANDOM_INTERLEAVING) && defined(ACC_LEVEL_INTERLEAVING)
            AU.addRequired<riAccCodeGen_ref::AccLevelRICodeGen_ref>();
#	elif defined(RANDOM_INTERLEAVING) && defined(ALL_LEVEL_INTERLEAVING)
			AU.addRequired<riCodeGen_ref::AllLevelRICodeGen_ref>();
#   elif defined(THREAD_TUNING)
            AU.addRequired<staticTuningCodeGen_ref::TuningCodeGen_ref>();
#   else 
            // AU.addRequired<modelCodeGen_ref::ModelCodeGen_ref>();
            AU.addRequired<modelCodeGenOpt_ref::ModelCodeGenOpt_ref>();
#   endif
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



