
add_llvm_loadable_module(LLVMSymFP
  symFP.cpp

  DEPENDS
  intrinsics_gen
  PLUGIN_TOOL
  opt
  )
