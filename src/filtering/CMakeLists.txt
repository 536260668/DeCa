include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/engine")
include_directories("${PROJECT_SOURCE_DIR}/src/tools")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")

aux_source_directory(. DIR_LIB_SRCS)
add_library(Mutect2FilteringEngineFunctions ${DIR_LIB_SRCS} ThresholdCalculator.cpp ThresholdCalculator.h FilterMutectCalls.cpp FilterMutectCalls.h FilteredHaplotypeFilter.cpp FilteredHaplotypeFilter.h BaseQualityFilter.cpp BaseQualityFilter.h MappingQualityFilter.cpp MappingQualityFilter.h)

target_link_libraries(Mutect2FilteringEngineFunctions ToolsFunctions)

