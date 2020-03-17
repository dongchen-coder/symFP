#ifndef ssCodeGen_ref_hpp
#define ssCodeGen_ref_hpp

#include <map>
#include <vector>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "sampleNumAnalysis.hpp"

using namespace llvm;
using namespace std;

#define CLS 64
#define DS 8

#define SAMPLING 2

// #define PARALLEL_CXX_THREAD

//#define PROFILE_SEARCH_REUSE

// #define SEARCH_REUSE_SAME_LOOP
//#define SEARCH_REUSE_DIFFERENT_LOOPS

#define DumpRTMR
#define REFERENCE_GROUP
// #define DumpRefLease

namespace ssCodeGen_ref {
    struct StaticSamplingCodeGen_ref : public FunctionPass {
        static char ID;
        StaticSamplingCodeGen_ref();

        map<Instruction*,  string> arrayName;
        map<Instruction*,  string> arrayExpression;

        uint64_t refGlobalNumber = 0;
        map<Instruction*, int> refNumber;

        map<Loop*, int> loopNumber;
        // map< string, int> refToSameArrayCnt;
        map<Value*,  string> indvName;

        map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, uint64_t> sampleNum;

        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initArrayName();

        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,  vector< string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);

        void rtHistoGen();
        void subBlkRTGen(); 
#ifdef DumpRTMR
        void rtDumpGen();
        void rtToMRGen();
        void mrDumpGen();
#elif defined(DumpRefLease)
        void accessRatioCalGen();
        void initHitsCostsGen();
        void getPPUCGen();
        void getMaxPPUCGen();
        void DumpRIGen();
        void RLGen();
#endif
#ifdef REFERENCE_GROUP
        void rtMergeGen();
#endif
        void headerGen();

        
        string getBound(Value* bound);
        string getBound_Start(Value* bound);
        
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> findLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops);
        
        void searchReuseDifferentLoopsUpdateFuncGen();
        void searchReuseDifferentLoopsCalFuncGen();
        bool searchReuseDifferentLoopsInitGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space);
        bool searchReuseDifferentLoopsBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space);
        
        void searchReuseSameLoopInitGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space);
        void searchReuseSameLoopBodyGen( string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, string space);
        
        bool refRTSearchGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree, bool GenFlag,  string refName, int useID,  vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> loops, vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> currentLoops, string space);
        void refRTGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void refRTBodyGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree,  string refName, int useID);
        
        void mainGen();
        
        bool runOnFunction(Function &F) override;
        void getAnalysisUsage(AnalysisUsage &AU) const override;
    };
}

#endif
