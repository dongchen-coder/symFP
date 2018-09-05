# Static Parallel Sampling (SPS)

This is the underdevelopment Static Parallel Sampling (SPS) repository (based on LLVM-4.0.0). 

SPS is aiming to provide detailed miss ratio curves statically. 
For each function, SPS consumes LLVM IR code as input. 
Then SPS analyzes the IR code and generates a specialized C++ static sampling code.
By compiling and running the generated C++ static sampling code, miss ratio curve will be derived.

License: Free for open source projects (GPL v.2.0).

Software artifact for PLDI 18 "Locality Analysis through Static Parallel Sampling" can be found in http://doi.org/10.5281/zenodo.1218771

## Usage

(1) Install LLVM 4.0.0

(2) git clone https://github.com/dongchen-coder/symFP.git to the path llvm-4.0.0.src/lib/Transforms/SymFP

(3) compile source to generate LLVMSymFP.dylib

(4) opt -load LLVMSymFP.dylib -sps -spsrate=0.01 \<TARGET.bc\> TARGET.bc.opt 2\> TARGET\_StaticSampling.cpp

(5) Compile TARGET\_StaticSampling.cpp with C++ compiler and run


## Code structure

There are 6 analysis passes, 2 code generation pass and 1 top wrapper pass.

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

### BranchAnalysis Pass
This pass is to extract the brach condition that associated with each reference

### GlobalVariableAnalysis Pass
This pass is to extract glabel variables

### ArgumentAnalysis Pass
This pass is to extract the argument information of a function

### IndexAnalysis Pass
This pass is to extract all the reference asscoiated with its index expression

### ssCodeGen Pass
This pass is to generate the static sampling code to derive static miss ratio curves. (Reference pair based approach)

### ssCodeGen_ref Pass
This pass is to generate the static sampling code to derive static miss ratio curves. (Reference based approach)


