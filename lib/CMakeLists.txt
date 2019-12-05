option(PLUM_ENABLE "open to enable PLUM to get sampler for OpenMP benchmarks" OFF)
add_library(LLVMSPS MODULE
  # List your source files here.
  # HBWMemOpAnalysis.hpp 
  # HBWMemOpAnalysis.cpp
  Utils.hpp
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

  loopTreeTransform.hpp
  loopTreeTransform.cpp

  samplingRate.hpp
  sampleNumAnalysis.hpp
  sampleNumAnalysis.cpp

  # ssCodeGen.hpp
  # ssCodeGen.cpp

  ssCodeGen_ref.hpp
  ssCodeGen_ref.cpp

  # rsCodeGen_ref.hpp
  # rsCodeGen_ref.cpp

  psCodeGen_ref.hpp
  psCodeGen_ref.cpp

  # riCodeGen_ref.hpp
  # riCodeGen_ref.cpp

  # uiCodeGen_ref.hpp
  # uiCodeGen_ref.cpp

  # metricCollectGen.hpp
  # metricCollectGen.cpp
)

# User C++11 to compile your pass (i.e., supply -std=c++11).
target_compile_features(LLVMSPS PRIVATE cxx_range_for cxx_auto_type)

# LLVM is (typically) built with no C++ RTTI. We need to match that;
# otherwise, we'll get linker errors about missing RTTI data.
set_target_properties(LLVMSPS PROPERTIES
   COMPILE_FLAGS "-fno-rtti"
)

# Comment this out to disable PLUM and switch it to SPS
add_definitions(-DPARALLEL)

# Define the interleaving pattern of PLUM. 
# All these macro will work if and only if its in PLUM mode (-DPARALLEL is defined)
# 
# - DRANDOM_INTERLEAVING:
#     Random Interleaving
# - DUNIFORM_INTERLEAVING:
#     1:1 Interleaving in iteration granularity
# By default is 1:1 interleaving in data granularity
# add_definitions(-DRANDOM_INTERLEAVING)
# add_definitions(-DUNIFORM_INTERLEAVING)


