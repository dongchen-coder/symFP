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

        struct ThreadNode: public TreeNodeBase {
            int threadNum;
        };

        /* LoopRefTree after the parallel transofrm */
        loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* PTLoopRefTree;

        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*>outMostLoops;

        void findAllOutMostLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);

        void insertThreadNode(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void tranverseLoopRefTree(TreeNodeBase* node);
#if defined(UNIFORM_INTERLEAVING)
        void UniformLoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
#elif defined(RANDOM_INTERLEAVING)
        void RandomLoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
#else
        void LoopTreeTransform(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
#endif
        void DumpLoopTree(TreeNodeBase* LTroot, string prefix);
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
