# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.20.5/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.20.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bigdreamerxixi/Mutect2Cpp-master

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug

# Include any dependencies generated for this target.
include src/utils/CMakeFiles/UtilsFunction.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/utils/CMakeFiles/UtilsFunction.dir/compiler_depend.make

# Include the progress variables for this target.
include src/utils/CMakeFiles/UtilsFunction.dir/progress.make

# Include the compile flags for this target's objects.
include src/utils/CMakeFiles/UtilsFunction.dir/flags.make

src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o: src/utils/CMakeFiles/UtilsFunction.dir/flags.make
src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o: ../src/utils/BaseUtils.cpp
src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o: src/utils/CMakeFiles/UtilsFunction.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o -MF CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o.d -o CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/utils/BaseUtils.cpp

src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/utils/BaseUtils.cpp > CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.i

src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/utils/BaseUtils.cpp -o CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.s

# Object files for target UtilsFunction
UtilsFunction_OBJECTS = \
"CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o"

# External object files for target UtilsFunction
UtilsFunction_EXTERNAL_OBJECTS =

src/utils/libUtilsFunction.a: src/utils/CMakeFiles/UtilsFunction.dir/BaseUtils.cpp.o
src/utils/libUtilsFunction.a: src/utils/CMakeFiles/UtilsFunction.dir/build.make
src/utils/libUtilsFunction.a: src/utils/CMakeFiles/UtilsFunction.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libUtilsFunction.a"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/UtilsFunction.dir/cmake_clean_target.cmake
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/UtilsFunction.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/utils/CMakeFiles/UtilsFunction.dir/build: src/utils/libUtilsFunction.a
.PHONY : src/utils/CMakeFiles/UtilsFunction.dir/build

src/utils/CMakeFiles/UtilsFunction.dir/clean:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils && $(CMAKE_COMMAND) -P CMakeFiles/UtilsFunction.dir/cmake_clean.cmake
.PHONY : src/utils/CMakeFiles/UtilsFunction.dir/clean

src/utils/CMakeFiles/UtilsFunction.dir/depend:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bigdreamerxixi/Mutect2Cpp-master /Users/bigdreamerxixi/Mutect2Cpp-master/src/utils /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/utils/CMakeFiles/UtilsFunction.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/utils/CMakeFiles/UtilsFunction.dir/depend

