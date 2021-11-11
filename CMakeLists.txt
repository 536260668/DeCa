cmake_minimum_required(VERSION 3.20)
project(Mutect2Cpp_master)

set(CMAKE_CXX_STANDARD 17)

# counterpart of gcc -I
include_directories(../htslib/htslib)
include_directories(src)

add_executable(Mutect2Cpp_master src/main.cpp src/SimpleInterval.cpp src/Mutect2Utils.cpp
                src/Locatable.cpp src/IntervalUtils.cpp src/ActivityProfileState.cpp src/AssemblyRegion.cpp
                src/AssemblyRegion.h src/SAMRecord.cpp src/SAMRecord.h src/main.cpp src/Mutect2Engine.cpp
                src/Mutect2Engine.h src/QualityUtils.cpp src/ReadFilter.cpp src/ReadFilter.h src/AlignmentContext.cpp
                src/AlignmentContext.h src/M2ArgumentCollection.h src/ReferenceCache.cpp src/ReferenceCache.h src/ActivityProfile.cpp
                src/ActivityProfile.h src/BandPassActivityProfile.cpp src/BandPassActivityProfile.h src/MathUtils.cpp src/MathUtils.h
                src/ReadPileup.cpp src/ReadPileup.h src/NaturalLogUtils.h src/NaturalLogUtils.cpp src/cache/DigammaCache.h src/cache/DigammaCache.cpp
                src/cache/Log10FactorialCache.h src/cache/Log10FactorialCache.cpp src/cache/IntToDoubleFunctionCache.h src/cache/IntToDoubleFunctionCache.cpp)

# counterpart of gcc -L
target_link_libraries(Mutect2Cpp_master ${PROJECT_SOURCE_DIR}/../htslib/libhts.so)

