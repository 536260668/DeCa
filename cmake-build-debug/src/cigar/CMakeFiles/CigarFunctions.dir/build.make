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
include src/cigar/CMakeFiles/CigarFunctions.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/cigar/CMakeFiles/CigarFunctions.dir/compiler_depend.make

# Include the progress variables for this target.
include src/cigar/CMakeFiles/CigarFunctions.dir/progress.make

# Include the compile flags for this target's objects.
include src/cigar/CMakeFiles/CigarFunctions.dir/flags.make

src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/flags.make
src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o: ../src/cigar/Cigar.cpp
src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o -MF CMakeFiles/CigarFunctions.dir/Cigar.cpp.o.d -o CMakeFiles/CigarFunctions.dir/Cigar.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/Cigar.cpp

src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CigarFunctions.dir/Cigar.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/Cigar.cpp > CMakeFiles/CigarFunctions.dir/Cigar.cpp.i

src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CigarFunctions.dir/Cigar.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/Cigar.cpp -o CMakeFiles/CigarFunctions.dir/Cigar.cpp.s

src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/flags.make
src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o: ../src/cigar/CigarElement.cpp
src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o -MF CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o.d -o CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarElement.cpp

src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CigarFunctions.dir/CigarElement.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarElement.cpp > CMakeFiles/CigarFunctions.dir/CigarElement.cpp.i

src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CigarFunctions.dir/CigarElement.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarElement.cpp -o CMakeFiles/CigarFunctions.dir/CigarElement.cpp.s

src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/flags.make
src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o: ../src/cigar/CigarOperator.cpp
src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o -MF CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o.d -o CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarOperator.cpp

src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarOperator.cpp > CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.i

src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/CigarOperator.cpp -o CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.s

src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/flags.make
src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o: ../src/cigar/TextCigarCodec.cpp
src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o: src/cigar/CMakeFiles/CigarFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o -MF CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o.d -o CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/TextCigarCodec.cpp

src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/TextCigarCodec.cpp > CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.i

src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar/TextCigarCodec.cpp -o CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.s

# Object files for target CigarFunctions
CigarFunctions_OBJECTS = \
"CMakeFiles/CigarFunctions.dir/Cigar.cpp.o" \
"CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o" \
"CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o" \
"CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o"

# External object files for target CigarFunctions
CigarFunctions_EXTERNAL_OBJECTS =

src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/Cigar.cpp.o
src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/CigarElement.cpp.o
src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/CigarOperator.cpp.o
src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/TextCigarCodec.cpp.o
src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/build.make
src/cigar/libCigarFunctions.a: src/cigar/CMakeFiles/CigarFunctions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libCigarFunctions.a"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && $(CMAKE_COMMAND) -P CMakeFiles/CigarFunctions.dir/cmake_clean_target.cmake
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CigarFunctions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/cigar/CMakeFiles/CigarFunctions.dir/build: src/cigar/libCigarFunctions.a
.PHONY : src/cigar/CMakeFiles/CigarFunctions.dir/build

src/cigar/CMakeFiles/CigarFunctions.dir/clean:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar && $(CMAKE_COMMAND) -P CMakeFiles/CigarFunctions.dir/cmake_clean.cmake
.PHONY : src/cigar/CMakeFiles/CigarFunctions.dir/clean

src/cigar/CMakeFiles/CigarFunctions.dir/depend:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bigdreamerxixi/Mutect2Cpp-master /Users/bigdreamerxixi/Mutect2Cpp-master/src/cigar /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar /Users/bigdreamerxixi/Mutect2Cpp-master/cmake-build-debug/src/cigar/CMakeFiles/CigarFunctions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/cigar/CMakeFiles/CigarFunctions.dir/depend

