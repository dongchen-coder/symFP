#ifndef loopTreeTransform_hpp
#define loopTreeTransform_hpp

#include <map>
#include <unordered_map>
#include <unordered_set>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "Utils.hpp"
#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "IVDependenceAnalysis.hpp"
#include "AccessGraphAnalysis.hpp"
#include "sampleNumAnalysis.hpp"

using namespace llvm;

typedef loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode LoopRefTNode;

namespace loopTreeTransform {
    struct ParallelLoopTreeTransform : public FunctionPass {
        static char ID;
        ParallelLoopTreeTransform();

        uint64_t total_cache_access;
        /* LoopRefTree after the parallel transofrm */
        LoopRefTNode* PTLoopRefTree;

        std::vector<LoopRefTNode*>outMostLoops;

        std::map<Instruction*, uint64_t> outMostLoopPerIterationSpace;
		
        /* all array references enclosed by each outermost looop */
        std::map<LoopRefTNode*, vector<Instruction*>> refPerOutMostLoop;

        void computePerIterationSpace();
        /* compute loop iteration space */
        uint64_t computeIterSpace(LoopRefTNode* LoopRefTree);

        void findAllOutMostLoops(LoopRefTNode* LoopRefTree);

        void tranverseLoopRefTree(LoopRefTNode* node);

        void LoopTreeTransform(LoopRefTNode* LoopRefTree);
        
        void DumpLoopTree(LoopRefTNode* LTroot, std::string prefix);
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
