include_directories("${PROJECT_SOURCE_DIR}/src")

aux_source_directory(. DIR_LIB_SRCS)

add_library(ReadFunctions ${DIR_LIB_SRCS} CigarUtils.cpp CigarUtils.h)