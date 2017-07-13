
add_llvm_loadable_module(LLVMSymFP
  symFP.cpp
  argAnalysis.cpp
  brchAnalysis.cpp
  gVarAnalysis.cpp
  idxAnalysis.cpp
  loopAnalysis.cpp
  

  DEPENDS
  intrinsics_gen
  PLUGIN_TOOL
  opt
  )
