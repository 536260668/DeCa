include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/engine")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")

aux_source_directory(. DIR_LIB_SRCS)
add_library(Mutect2FilteringEngineFunctions ${DIR_LIB_SRCS} M2FiltersArgumentCollection.cpp M2FiltersArgumentCollection.h StrandArtifactFilter.cpp StrandArtifactFilter.h ErrorProbabilities.cpp ErrorProbabilities.h)

