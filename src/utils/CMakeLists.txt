aux_source_directory(. DIR_LIB_SRCS)
add_library(UtilsFunctions ${DIR_LIB_SRCS} PeUtils.cpp PeUtils.h LocusIteratorByState.cpp LocusIteratorByState.h)