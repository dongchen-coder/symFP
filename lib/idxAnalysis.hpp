//
//  idxAnalysis.hpp
//  LLVM
//
//  Created by dongchen on 5/18/17.
//
//



#ifndef idxAnalysis_hpp
#define idxAnalysis_hpp

#include <stdio.h>
#include "llvm/Pass.h"
#include "llvm/IR/Function.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/IR/InstIterator.h"
#include "llvm/IR/Instructions.h"
#include <map>
#include <vector>

/* Debugging flag, enable if debug is needed */
//#define IDX_DEBUG

using namespace llvm;

/* This is a function pass that extracts array accesses information (reference name and its indexing function) */

namespace idxAnalysis {
    struct IndexAnalysis : public FunctionPass {
        static char ID;
        IndexAnalysis();
        
        /* Using hash table to associating the load/store instructions with the array names and expression */
        std::map<Instruction*, std::string> arrayName;
        std::map<Instruction*, std::string> arrayExpression;
        std::map<Instruction*, std::vector<std::string>> arrayAccessVariable;

        /* 
         * Using hash table to associating the array load/store instruction with its cost value 
         * now the table is reference based, 2 accesses to the same array will compute its cost separately
         * Will merge at the end of this pass -> BCArrayList;
         * assign the Expected BCArray ratio to a global value -> BCArrayRatio
         */
        std::map<Instruction*, int> BCReference;
        std::map<std::string, int> BCArrayList;

        /* Expeted cost */
        double BCArrayRatio;
        
        /* Extracting array index expression from load/store instruction */
        std::vector<Instruction*> instStackInit(Instruction* inst);
        std::string computeExpression(Instruction* inst, std::vector<std::string>& idvs);
        
        /* Extracting array name from load/store instruction */
        std::string getArrayName(GetElementPtrInst* inst);
        
        /* Scan all IR instructions in function to find load/stores accessing arrays */
        void findAllArrayAccesses(Function &F);
        
        /* Dump (array name, index expression) pairs for all references inside the function */
        void dumpAllInfo();
        
        /* Analysis pass main function */
        bool runOnFunction(Function &F) override;
    };
}


#endif /* idxAnalysis_hpp */
