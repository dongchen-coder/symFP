#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"

using namespace llvm;

#define RANDOM_REF_SAMPLING_RATE 0.02

namespace sampleNumAnalysis {
    struct SampleNumberAnalysis : public FunctionPass {
        static char ID;
        SampleNumberAnalysis();
        
        std::map<Instruction*, int> refNumber;
        std::map<Loop*, int> loopNumber;
        std::map<Value*, std::string> indvName;
        
        int getConstantStride(Value* expr);
        
        string getBound(Value* bound);
        
        std::map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, uint64_t> sampleNum;
        
        int valueOfExpression(Value* expr, std::map<Value*, int> counter);
        
        void calculateSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        
        void dumpSampleNum(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LTroot, std::string prefix);
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
        
    };
}
