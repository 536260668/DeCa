cmake_minimum_required(VERSION 3.20)
project(Mutect2Cpp_master)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_BUILD_TYPE "Release")
#
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")

#set(CMAKE_C_FLAGS_Release "${CMAKE_C_FLAGS_Release} -O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -O3 -fno-omit-frame-pointer -g")

# counterpart of gcc -I
include_directories(../htslib)
include_directories(src)

aux_source_directory(./src DIR_SRCS)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/src/cache")
include_directories("${PROJECT_SOURCE_DIR}/src/cigar")
include_directories("${PROJECT_SOURCE_DIR}/src/haplotype")
include_directories("${PROJECT_SOURCE_DIR}/src/path")
include_directories("${PROJECT_SOURCE_DIR}/src/read")
include_directories("${PROJECT_SOURCE_DIR}/src/graph")
include_directories("${PROJECT_SOURCE_DIR}/src/intel/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/smithwaterman")
include_directories("${PROJECT_SOURCE_DIR}/src/Mutect2")
include_directories("${PROJECT_SOURCE_DIR}/src/variantcontext")

add_subdirectory(./src/cache)
add_subdirectory(./src/cigar)
add_subdirectory(./src/haplotype)
add_subdirectory(./src/path)
add_subdirectory(./src/read)
add_subdirectory(./src/graph)
add_subdirectory(./src/intel/smithwaterman)
add_subdirectory(./src/smithwaterman)
add_subdirectory(./src/Mutect2)
add_subdirectory(./src/variantcontext)
add_subdirectory(./src/haplotypecaller)
add_subdirectory(./src/utils)
add_subdirectory(./src/tools)
add_subdirectory(./src/samtools)
add_subdirectory(./src/param)
add_subdirectory(./src/clipping)

add_executable(Mutect2Cpp_master ${DIR_SRCS} src/engine/ReferenceContext.cpp src/engine/ReferenceContext.h src/engine/AlignmentContext.cpp src/engine/AlignmentContext.h src/transfer/PalindromeArtifactClipReadTransformer.cpp src/transfer/PalindromeArtifactClipReadTransformer.h src/pileRead.h src/collection/IndexedSet.h src/collection/IndexedAlleleList.h)


# counterpart of gcc -L
target_link_libraries(Mutect2Cpp_master ${PROJECT_SOURCE_DIR}/../htslib/libhts.dylib)
target_link_libraries(Mutect2Cpp_master CigarFunctions)
target_link_libraries(Mutect2Cpp_master ReadFunctions)
target_link_libraries(Mutect2Cpp_master CacheFunctions)
target_link_libraries(Mutect2Cpp_master GraphFunctions)
target_link_libraries(Mutect2Cpp_master HaplotypeFunctions)
target_link_libraries(Mutect2Cpp_master PathFunctions)
target_link_libraries(Mutect2Cpp_master smithwatermanFunctions)
target_link_libraries(Mutect2Cpp_master SWFunctions)
target_link_libraries(Mutect2Cpp_master Mutect2Functions)
target_link_libraries(Mutect2Cpp_master VariantContext)
target_link_libraries(Mutect2Cpp_master HaplotypeCallerFunctions)
target_link_libraries(Mutect2Cpp_master UtilsFunctions)
target_link_libraries(Mutect2Cpp_master ToolsFunctions)
target_link_libraries(Mutect2Cpp_master SamtoolsFunctions)
target_link_libraries(Mutect2Cpp_master ParamFunctions)
target_link_libraries(Mutect2Cpp_master ClippingFunctions)



