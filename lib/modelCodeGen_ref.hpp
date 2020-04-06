#ifndef modelCodeGen_ref_hpp
#define modelCodeGen_ref_hpp

#include <map>
#include <vector>

#include "llvm/Pass.h"
#include "llvm/IR/Function.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "loopTreeTransform.hpp"
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
// #define REFERENCE_GROUP
// #define DumpRefLease

namespace modelCodeGen_ref {
    enum RefType
    {
        OUTERMOST_ONLY, // A[I]
        OUTERMOST_INVARIANT, // A[J]
        REVERSE, // A[J][I]
        REGULAR // A[I][J]
    };
    struct ModelCodeGen_ref : public FunctionPass {
        static char ID;
        ModelCodeGen_ref();

        map<Instruction*,  string> arrayName;
        map<Instruction*,  string> arrayExpression;
        map<Instruction*, std::vector<std::string>> arrayAccessVariable;
        map<string, RefType>arrayTypeMap;
        map<string, uint64_t> outMostLoopPerIterationSpace;

        uint64_t refGlobalNumber = 0;
        map<Instruction*, int> refNumber;

        map<Loop*, int> loopNumber;
        // map< string, int> refToSameArrayCnt;
        map<Value*,  string> indvName;

        map<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*, uint64_t> sampleNum;
        vector<loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*> outloops;

        void numberRefToSameArray(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void numberLoops(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);
        void initIndvName(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        void initArrayName();

        void defineRefs(Instruction* arrayInstr);

        void addrCalFuncGen(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree,  vector< string> indvs);
        void addrCalFuncGenTop(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode *LoopRefTree);

        void rtHistoGen();
        void subBlkRTGen(); 

        void parallelModelCodeGen();
#ifdef DumpRTMR
        void rtDumpGen();
        void rtToMRGen();
        void mrDumpGen();
#endif

#ifdef REFERENCE_GROUP
        void rtMergeGen();
#endif
        void headerGen();

        uint64_t computeIterSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode* LoopRefTree);
        uint64_t getOuterMostLoopIterationSpace(loopAnalysis::LoopIndvBoundAnalysis::LoopRefTNode*);
        string getBound(Value* bound);
        string getBound_Start(Value* bound);
        string getLoopInc(Value *inc);
        
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
