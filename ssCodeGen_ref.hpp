#ifndef ssCodeGen_ref_hpp
#define ssCodeGen_ref_hpp

#include <map>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"

using namespace llvm;

#define CLS 32
#define DS 8

#define SAMPLING 2
#define RANDOM_REF_SAMPLING_RATE 0.02

namespace ssCodeGen_ref {
    struct StaticSamplingCodeGen_ref : public FunctionPass {
        static char ID;
        StaticSamplingCodeGen_ref();
        
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        std::map<Instruction*, int> refNumber;
        std::map<Loop*, int> loopNumber;
        std::map<std::string, int> refToSameArrayCnt;
        std::map<Value*, std::string> indvName;
        
        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initArrayName();
        
        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree, std::vector<std::string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        
        void rtHistoGen();
        void rtDumpGen();
        void rtToMRGen();
        void mrDumpGen();
        void headerGen();
        
        string getBound(Value *bound);
        std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        
        bool refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag, std::string refName, int useID, std::vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space);
        void refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, std::string refName, int useID);
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
