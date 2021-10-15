cmake_minimum_required(VERSION 3.20)
project(Mutect2Cpp_master)

set(CMAKE_CXX_STANDARD 14)

# counterpart of gcc -I
include_directories(../htslib/htslib)
include_directories(src)

<<<<<<< Updated upstream
add_executable(Mutect2Cpp_master src/main.cpp src/SimpleInterval.cpp src/Mutect2Utils.cpp src/Locatable.cpp src/IntervalUtils.cpp src/SAMRecord.cpp src/SAMRecord.h)

# counterpart of gcc -L
target_link_libraries(Mutect2Cpp_master /home/lhh/Documents/learning/project/Mutect2/htslib/libhts.so)
=======
add_executable(Mutect2Cpp_master src/main.cpp src/SimpleInterval.cpp src/Mutect2Utils.cpp src/Locatable.cpp src/IntervalUtils.cpp src/ActivityProfileState.cpp src/AssemblyRegion.cpp src/AssemblyRegion.h)
#add_executable(Mutect2Cpp_master src/main.cpp src/ActivityProfileState.cpp)
>>>>>>> Stashed changes
