include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/Mutect2FilteringEngine")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")
include_directories("${PROJECT_SOURCE_DIR}/src/Mutect2")

aux_source_directory(. DIR_LIB_SRCS)
add_library(clusteringFunctions ${DIR_LIB_SRCS} SequencingError.cpp SequencingError.h BetaBinomialCluster.cpp BetaBinomialCluster.h BetaDistributionShape.cpp BetaDistributionShape.h)