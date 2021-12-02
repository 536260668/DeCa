include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/haplotype")
include_directories("${PROJECT_SOURCE_DIR}/src/path")
include_directories("${PROJECT_SOURCE_DIR}/src/read")
include_directories("${PROJECT_SOURCE_DIR}/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/src/intel/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/smithwaterman")

aux_source_directory(. DIR_LIB_SRCS)

add_library(Mutect2Functions ${DIR_LIB_SRCS} AssemblyResultSet.cpp AssemblyResultSet.h)