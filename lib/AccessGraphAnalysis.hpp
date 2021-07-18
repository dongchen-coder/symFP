//
//  AccessGraphAnalysis.hpp
//  LLVMSPS
//
//  Created by noya-fangzhou on 5/10/21.
//

#ifndef AccessGraphAnalysis_hpp
#define AccessGraphAnalysis_hpp

#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <string>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Analysis/LoopInfo.h"
#include "llvm/Analysis/CFG.h"
#include "llvm/IR/Constants.h"

#include "idxAnalysis.hpp"
#include "argAnalysis.hpp"
#include "gVarAnalysis.hpp"
#include "loopAnalysis.hpp"
#include "plumCodeGenUtil.hpp"

using namespace llvm;
using namespace std;

namespace AccGraphAnalysis {
	struct AccessGraphAnalysis : public FunctionPass {
		static char ID;
		AccessGraphAnalysis();
		
		unordered_map<LoopRefTNode *, AccGraph *> AccGraphLoopMap;
		
		void buildAccessGraphUtil(LoopRefTNode * L);
		void buildAccessGraph(LoopRefTNode * L);
		void dumpAccessGraph(AccGraph * G, string prefix);
		
		/* Analysis pass main function */
		bool runOnFunction(Function &F) override;

		void getAnalysisUsage(AnalysisUsage &AU) const override;
	};
}
#endif /* AccessGraphAnalysis_hpp */
