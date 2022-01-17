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
CMAKE_BINARY_DIR = /Users/bigdreamerxixi/Mutect2Cpp-master/build

# Include any dependencies generated for this target.
include src/haplotype/CMakeFiles/HaplotypeFunctions.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/haplotype/CMakeFiles/HaplotypeFunctions.dir/compiler_depend.make

# Include the progress variables for this target.
include src/haplotype/CMakeFiles/HaplotypeFunctions.dir/progress.make

# Include the compile flags for this target's objects.
include src/haplotype/CMakeFiles/HaplotypeFunctions.dir/flags.make

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/flags.make
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o: ../src/haplotype/Allele.cpp
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o -MF CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o.d -o CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Allele.cpp

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Allele.cpp > CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.i

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Allele.cpp -o CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.s

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/flags.make
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o: ../src/haplotype/EventMap.cpp
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o -MF CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o.d -o CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/EventMap.cpp

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/EventMap.cpp > CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.i

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/EventMap.cpp -o CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.s

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/flags.make
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o: ../src/haplotype/Haplotype.cpp
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o -MF CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o.d -o CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Haplotype.cpp

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Haplotype.cpp > CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.i

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype/Haplotype.cpp -o CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.s

# Object files for target HaplotypeFunctions
HaplotypeFunctions_OBJECTS = \
"CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o" \
"CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o" \
"CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o"

# External object files for target HaplotypeFunctions
HaplotypeFunctions_EXTERNAL_OBJECTS =

src/haplotype/libHaplotypeFunctions.a: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Allele.cpp.o
src/haplotype/libHaplotypeFunctions.a: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/EventMap.cpp.o
src/haplotype/libHaplotypeFunctions.a: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/Haplotype.cpp.o
src/haplotype/libHaplotypeFunctions.a: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/build.make
src/haplotype/libHaplotypeFunctions.a: src/haplotype/CMakeFiles/HaplotypeFunctions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libHaplotypeFunctions.a"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && $(CMAKE_COMMAND) -P CMakeFiles/HaplotypeFunctions.dir/cmake_clean_target.cmake
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/HaplotypeFunctions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/haplotype/CMakeFiles/HaplotypeFunctions.dir/build: src/haplotype/libHaplotypeFunctions.a
.PHONY : src/haplotype/CMakeFiles/HaplotypeFunctions.dir/build

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/clean:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype && $(CMAKE_COMMAND) -P CMakeFiles/HaplotypeFunctions.dir/cmake_clean.cmake
.PHONY : src/haplotype/CMakeFiles/HaplotypeFunctions.dir/clean

src/haplotype/CMakeFiles/HaplotypeFunctions.dir/depend:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bigdreamerxixi/Mutect2Cpp-master /Users/bigdreamerxixi/Mutect2Cpp-master/src/haplotype /Users/bigdreamerxixi/Mutect2Cpp-master/build /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/haplotype/CMakeFiles/HaplotypeFunctions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/haplotype/CMakeFiles/HaplotypeFunctions.dir/depend

