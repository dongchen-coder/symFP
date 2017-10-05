# Symbolic Footprint

This is the underdevelopment symblic footprint repository (based on LLVM-4.0.0). 

Symbolic footprint is aiming to provide detailed miss ratio curves statically for each function.

## Code structure

There are 5 analysis passes, 1 code generation pass and 1 top wrapper pass. 

### symFP Pass (Top)

symFP is the top wrapper analysis pass which runs other 5 analysis passes (ArgumentAnalysis, BranchAnalysis, GlobalVariableAnalysis, IndexAnalysis, LoopIndvBoundAnalysis) and 1 code generation pass (ssCodeGen). 

### LoopIndvBoundAnalysis Pass 

This analysis abstract the function as a tree contains non-leaf loop node (AA = NULL) and leaf reference node (L = NULL) by LoopRefTNode.   

```C++
struct LoopRefTNode {
            int LoopLevel;
            Loop* L;
            LoopInfoStruct* LIS;
            Instruction* AA;
            vector<LoopRefTNode *>* next;
};
```

### BranchAnalysis Pass (Not finished)
This pass is to extract the brach condition that associated with each reference

### GlobalVariableAnalysis Pass
This pass is to extract glabel variables

### ArgumentAnalysis Pass
This pass is to extract the argument information of a function

### IndexAnalysis Pass
This pass is to extract all the reference asscoiated with its index expression

### ssCodeGen Pass
This pass is to generate the static sampling code to derive static miss ratio curves.


