#ifndef loopTreeTransform_hpp
#define loopTreeTransform_hpp

#include <map>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "Utils.hpp"
#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "sampleNumAnalysis.hpp"

using namespace llvm;

namespace loopTreeTransform {
    struct ParallelLoopTreeTransform : public FunctionPass {
        static char ID;
        ParallelLoopTreeTransform();

        typedef loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode LoopRefTNode;

        /* LoopRefTree after the parallel transofrm */
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* PTLoopRefTree;

        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>outMostLoops;

        std::map<Instruction*, uint64_t> outMostLoopPerIterationSpace;


        void computePerIterationSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);

        void findAllOutMostLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);

        void tranverseLoopRefTree(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* node);

        void LoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        
        void DumpLoopTree(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LTroot, string prefix);
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
