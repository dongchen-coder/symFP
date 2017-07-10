# Symbolic Footprint

This is the underdevelopment symblic footprint repository. 

symFP is the major analysis pass which is dependent on the other four analysis passes (index, argument, global variable and loop information analysis pass).

## LoopAnalysis Pass
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
### Induction Variable
Iterate each Condition BasicBlocks in each Loop and find the <i>ICMPInstr</i>, then the first operand is its induction variable

### Loop Bound
Iterating each Condition BasicBlocks in each Loop and find the <i>ICMPInstr</i>, then the second operand is its upper bound. Meanwhile, iterating the predecessors of this Condition BasicBlock, ignoring the Incremental BasicBlocks, then go over the <i>STOREInstr</i>, when the second operand of the <i>STOREInstr</i> is equal to the current Induction Variable, then we say the first operand is the lower bound of the current Loop.




