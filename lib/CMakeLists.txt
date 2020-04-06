option(PLUM_ENABLE "open to enable PLUM to get sampler for OpenMP benchmarks" OFF)
add_library(LLVMSPS MODULE
  # List your source files here.
  # HBWMemOpAnalysis.hpp 
  # HBWMemOpAnalysis.cpp
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

  modelCodeGen_ref.hpp
  modelCodeGen_ref.cpp

  # rsCodeGen_ref.hpp
  # rsCodeGen_ref.cpp

  riIterCodeGen_ref.hpp
  riIterCodeGen_ref.cpp

  riAccCodeGen_ref.hpp
  riAccCodeGen_ref.cpp

  # riCodeGen_ref.hpp
  # riCodeGen_ref.cpp

  uiIterCodeGen_ref.hpp
  uiIterCodeGen_ref.cpp

  uiAccCodeGen_ref.hpp
  uiAccCodeGen_ref.cpp

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
# set_target_properties(LLVMSPS PROPERTIES
#    LINK_FLAGS "-undefined dynamic_lookup"
# )


# Comment this out to disable PLUM and switch it to SPS
add_definitions(-DPARALLEL)

# Define the interleaving pattern of PLUM. 
# All these macro will work if and only if its in PLUM mode (-DPARALLEL is defined)
# 
# - DRANDOM_INTERLEAVING:
#     Random Interleaving
# - DUNIFORM_INTERLEAVING:
#     Uniform Interleaving
# By default is uniform interleaving
# add_definitions(-DRANDOM_INTERLEAVING)
add_definitions(-DUNIFORM_INTERLEAVING)

# - DACC_LEVEL_INTERLEAVING:
#     Do Interleaving per each access
# - DITER_LEVEL_INTERLEAVING:
#     Do Interleaving per each iteration
# - DALL_LEVEL_INTERLEAVING:
#     Do interleaving for all references and interations
# By default is 1:1 interleaving in data granularity
# add_definitions(-DITER_LEVEL_INTERLEAVING)
add_definitions(-DACC_LEVEL_INTERLEAVING)
# add_definitions(-DALL_LEVEL_INTERLEAVING)

# - DUNIFORM_SMOOTHING
#     Do the RI calibration in uniform distribution
# - DGAUSSIAN_SMOOTHING #     Do the RI calibration in gaussian distribution
# By default is uniform smoothing
# add_definitions(-DGAUSSIAN_SMOOTHING)
# add_definitions(-DUNIFORM_SMOOTHING)



