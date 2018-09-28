# Static Parallel Sampling (SPS)

This is the **Underdevelopment** Static Parallel Sampling (SPS) repository (based on LLVM-6.0.0). We keep on improving this tool and introducing new techniques in this repository.

Software artifact for PLDI 18 "Locality Analysis through Static Parallel Sampling" can be found in http://doi.org/10.5281/zenodo.1218771

SPS is aiming to provide detailed miss ratio curves statically and efficiently. 
For each function, SPS consumes LLVM IR code as input. 
Then SPS analyzes the IR code and generates a specialized C++ static sampling code.
By compiling and running the generated C++ static sampling code, miss ratio curve will be derived.

License: Free for open source projects (GPL v.2.0).

## Usage

(1) Install LLVM 6.0.0

(2) git clone https://github.com/dongchen-coder/symFP.git to the path llvm-6.0.0.src/lib/Transforms/SymFP

(3) compile source to generate LLVMStaticSampling.dylib

(4) opt -load LLVMStaticSampling.dylib -sps -spsrate=0.01 \<TARGET.bc\> TARGET.bc.opt 2\> TARGET\_StaticSampling.cpp

(5) Compile TARGET\_StaticSampling.cpp with C++ compiler and run

NOTE: we also provides CMake files for standalone build for linux named (CmakeLists.txt.standalone).

IMPORTANT: our tool currently relys on DEBUG build of LLVM libs and opt tool. We recommand you build llvm from source with DEBUG mode before you start to play with our tool.

## Code structure

There are 6 analysis passes, 2 code generation pass and 1 top wrapper pass.

### Static Sampling Pass (sps.cpp)

Static sampling pass is the top wrapper analysis pass which runs other 5 analysis passes (ArgumentAnalysis, BranchAnalysis, GlobalVariableAnalysis, IndexAnalysis, LoopIndvBoundAnalysis) and 1 code generation pass (ssCodeGen). 

### Argument Analysis Pass (argAnalysis.cpp)

Argument analysis pass is to extract the argument information of a function

### Branch Analysis Pass (brchAnalysis.cpp, disabled, need more work)

Branch analysis pass is to extract the brach condition that associated with each reference

### Global Variable Analysis Pass (gVarAnalysis.cpp)

Global variable analysis pass is to extract glabel variables

### Index Analysis Pass (idxAnalysis.cpp)

Index analysis pass is to extract all the reference asscoiated with its index expression

### Loop Analysis Pass (loopAnalysis.cpp)

Loop analysis pass abstracts the function as a tree contains non-leaf loop node (AA = NULL) and leaf reference node (L = NULL) by LoopRefTNode.   

```C++
struct LoopRefTNode {
            int LoopLevel;
            Loop* L;
            LoopInfoStruct* LIS;
            Instruction* AA;
            vector<LoopRefTNode *>* next;
};
```
### Sample Number Analysis Pass (sampleNumAnalysis.cpp)

Sample number analysis pass calculates the number of samples needed according to sampling rate for each loop.

### Static Sampling Code Generating Pass (ssCodeGen.cpp, reference pair based)

This pass is to generate the static sampling code to derive static miss ratio curves. (Reference pair based approach)

### Static Sampling Code Generating  Pass (ssCodeGen.cpp, reference based)

This pass is to generate the static sampling code to derive static miss ratio curves. (Reference based approach)


