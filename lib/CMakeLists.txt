option(ENABLE_PLUM
    "open to enable PLUM to get sampler for OpenMP benchmarks" 
    OFF
)
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

  modelCodeGenOpt_ref.hpp
  modelCodeGenOpt_ref.cpp
  # rsCodeGen_ref.hpp
  # rsCodeGen_ref.cpp

  riIterCodeGen_ref.hpp
  riIterCodeGen_ref.cpp

  riAccCodeGen_ref.hpp
  riAccCodeGen_ref.cpp

  riCodeGen_ref.hpp
  riCodeGen_ref.cpp

  uiIterCodeGen_ref.hpp
  uiIterCodeGen_ref.cpp

  uiAccCodeGen_ref.hpp
  uiAccCodeGen_ref.cpp
 
  uiAccCodeGenOpt_ref.hpp
  uiAccCodeGenOpt_ref.cpp
  
  staticTuningCodeGen_ref.hpp
  staticTuningCodeGen_ref.cpp

  IVDependenceAnalysis.hpp
  IVDependenceAnalysis.cpp

  dynSchCodeGen_ref.hpp
  dynSchCodeGen_ref.cpp

  statSchCodeGen_ref.hpp
  statSchCodeGen_ref.cpp

  plumCodeGen_ref.hpp
  plumCodeGen_ref.cpp

  plumSamplerCodeGenerator.hpp
  plumSamplerCodeGenerator.cpp
  
  plumCodeGenUtil.hpp
  plumCodeGenUtil.cpp

  AccessGraphAnalysis.hpp
  AccessGraphAnalysis.cpp
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


# Define the interleaving pattern of PLUM. 
# All these macro will work if and only if its in PLUM mode (-DPARALLEL is defined)
# 
# - DRANDOM_INTERLEAVING:
#     Random Interleaving
# - DUNIFORM_INTERLEAVING:
#     Uniform Interleaving
#
# - DSTATIC_SCHEDULING:
#     running in OpenMP static scheduling
# - DDYNAMIC_SCHEDULING:
#     running in OpenMP dynamic scheduling

OPTION(ENABLE_UNIFORM_INTERLEAVING
  "Use the Uniform Interleaving"
  OFF
)

OPTION(ENABLE_RANDOM_INTERLEAVING
  "Use the Random Interleaving"
  OFF
)

OPTION(STATIC_SCHEDULING
  "In OpenMP static scheduling"
  OFF
)

OPTION(DYNAMIC_SCHEDULING
  "In OpenMP dynamic scheduling"
  OFF
)

OPTION(ENABLE_PARALLEL_SAMPLING
  "Open the Paralle Sampling Optimization"
  OFF
)

if(ENABLE_PLUM)
    # Comment this out to disable PLUM and switch it to SPS
    add_definitions(-DPARALLEL)
endif(ENABLE_PLUM)

if(STATIC_SCHEDULING)
  add_definitions(-DSTATIC_SCHEDULING)
elseif(DYNAMIC_SCHEDULING) 
  add_definitions(-DDYNAMIC_SCHEDULING)
endif(STATIC_SCHEDULING)

if(ENABLE_UNIFORM_INTERLEAVING)
  add_definitions(-DUNIFORM_INTERLEAVING)
  add_definitions(-DACC_LEVEL_INTERLEAVING)
endif(ENABLE_UNIFORM_INTERLEAVING)


if(ENABLE_RANDOM_INTERLEAVING)
  add_definitions(-DRANDOM_INTERLEAVING)
  add_definitions(-DALL_LEVEL_INTERLEAVING)
endif(ENABLE_RANDOM_INTERLEAVING)


# enable parallel sampling optimization
if(ENABLE_PARALLEL_SAMPLING)
  add_definitions(-DPARALLEL_CXX_THREAD)
endif(ENABLE_PARALLEL_SAMPLING)
