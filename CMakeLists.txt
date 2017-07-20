
add_llvm_loadable_module(LLVMSymFP
  symFP.cpp
  argAnalysis.hpp
  argAnalysis.cpp

  brchAnalysis.hpp
  brchAnalysis.cpp

  gVarAnalysis.hpp
  gVarAnalysis.cpp

  idxAnalysis.hpp
  idxAnalysis.cpp

  loopAnalysis.hpp
  loopAnalysis.cpp
  

  DEPENDS
  intrinsics_gen
  PLUGIN_TOOL
  opt
  )
