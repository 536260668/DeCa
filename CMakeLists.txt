cmake_minimum_required(VERSION 3.20)
project(Mutect2Cpp_master)

set(CMAKE_CXX_STANDARD 14)

include_directories(src)

add_executable(Mutect2Cpp_master src/main.cpp src/SimpleInterval.cpp src/Mutect2Utils.cpp src/Locatable.cpp src/IntervalUtils.cpp)
