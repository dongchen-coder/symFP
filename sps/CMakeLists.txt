add_library(LLVMsps MODULE
  sps.cpp

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

  samplingRate.hpp
  sampleNumAnalysis.hpp
  sampleNumAnalysis.cpp

  ssCodeGen.hpp
  ssCodeGen.cpp

  ssCodeGen_ref.hpp
  ssCodeGen_ref.cpp
  )

target_compile_features(LLVMsps PRIVATE cxx_range_for cxx_auto_type)

set_target_properties(LLVMsps PROPERTIES
	COMPILE_FLAGS "-fno-rtti"
)

set_target_properties(LLVMsps PROPERTIES
	LINK_FLAGS "-undefined dynamic_lookup"
)

