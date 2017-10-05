# Symbolic Footprint

This is the underdevelopment symblic footprint repository (based on LLVM-4.0.0). 

Symbolic footprint is aiming to provide detailed miss ratio curves statically for each function.

## Code structure

There are 5 analysis passes, 1 code generation pass and 1 top wrapper pass. 

### symFP Pass (Top)

symFP is the top wrapper analysis pass which runs other 5 analysis passes (ArgumentAnalysis, BranchAnalysis, GlobalVariableAnalysis, IndexAnalysis, LoopIndvBoundAnalysis) and 1 code generation pass (ssCodeGen). 

### LoopIndvAnalysis Pass 

This pass analyze the loop information in the program. All Loop information was encapsulated in a LoopInfoStruct structrue.
Here is the detail:
```C++

 struct LoopInfoStruct{
        Loop *L;                              // Loop Object  
        Value * IDV;                          // Induction Var in the Loop
        LoopIndvBoundAnalysis::LoopBound LB;  // The Loop Bound
        vector<Loop *> SL;                    // The SubLoops it contains
 };
```
#### Induction Variable
Iterate each Condition BasicBlocks in each Loop and find the <i>ICMPInst</i>, then the first operand is its induction variable

#### Loop Bound
Iterating each Condition BasicBlocks in each Loop and find the <i>ICMPInst</i>, then the second operand is its upper bound. Meanwhile, iterating the predecessors of this Condition BasicBlock, ignoring the Incremental BasicBlocks, then go over the <i>STOREInst</i>, when the second operand of the <i>STOREInst</i> is equal to the current Induction Variable, then we say the first operand is the lower bound of the current Loop.

### LoopBranchAnalysis Pass
This pass analyze the loop information, which contains the if-statement in the loop body. All Loop Branch Information was encapsulated in a LoopBranchStruct structure.
Here is the detail:
```C++

 struct LoopBranchStruct{
        Loop *L;                              // Loop Object  
        vector<ICmpInst *> BR;                // The branch condition within the loop
        vector<Value *> ARR;                  // The array within the branch
 };
```
#### Branch Condition
Iterate each BasicBlock in each Loop and extract those BasicBlocks belongs to the if-body, which are located between the BasicBlock named if.then and if.end. The Branch Condition is hidden in the predecessor BasicBlock of the begining of the if-body.

#### Valid Branch Condition
When locating the Branch Condition, we should then iterate each Instructions within the if-body, checking whether there's any array visit (<i>GetElementPtrInst</i>) in the if-body. In addition, we should check whether the variable in the Branch Condition contains the loop induction variable, which was found in the LoopIndvAnalysis Pass.


