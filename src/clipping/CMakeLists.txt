aux_source_directory(. DIR_LIB_SRCS)
add_library(ClippingFunctions ${DIR_LIB_SRCS} ReadClipper.cpp ReadClipper.h ClippingOp.cpp ClippingOp.h)