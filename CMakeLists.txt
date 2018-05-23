
add_llvm_loadable_module(LLVMSymFP
  sps/symFP.cpp
  sps/argAnalysis.hpp
  sps/argAnalysis.cpp

  sps/brchAnalysis.hpp
  sps/brchAnalysis.cpp

  sps/gVarAnalysis.hpp
  sps/gVarAnalysis.cpp

  sps/idxAnalysis.hpp
  sps/idxAnalysis.cpp

  sps/loopAnalysis.hpp
  sps/loopAnalysis.cpp

  sps/sampleNumAnalysis.hpp
  sps/sampleNumAnalysis.cpp

  sps/ssCodeGen.hpp
  sps/ssCodeGen.cpp

  sps/ssCodeGen_ref.hpp
  sps/ssCodeGen_ref.cpp

  DEPENDS
  intrinsics_gen
  PLUGIN_TOOL
  opt
  )
