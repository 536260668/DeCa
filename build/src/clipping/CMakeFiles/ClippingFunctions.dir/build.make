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
include src/clipping/CMakeFiles/ClippingFunctions.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/clipping/CMakeFiles/ClippingFunctions.dir/compiler_depend.make

# Include the progress variables for this target.
include src/clipping/CMakeFiles/ClippingFunctions.dir/progress.make

# Include the compile flags for this target's objects.
include src/clipping/CMakeFiles/ClippingFunctions.dir/flags.make

src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o: src/clipping/CMakeFiles/ClippingFunctions.dir/flags.make
src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o: ../src/clipping/ClippingOp.cpp
src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o: src/clipping/CMakeFiles/ClippingFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o -MF CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o.d -o CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ClippingOp.cpp

src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ClippingOp.cpp > CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.i

src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ClippingOp.cpp -o CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.s

src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o: src/clipping/CMakeFiles/ClippingFunctions.dir/flags.make
src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o: ../src/clipping/ReadClipper.cpp
src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o: src/clipping/CMakeFiles/ClippingFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o -MF CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o.d -o CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ReadClipper.cpp

src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ReadClipper.cpp > CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.i

src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping/ReadClipper.cpp -o CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.s

# Object files for target ClippingFunctions
ClippingFunctions_OBJECTS = \
"CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o" \
"CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o"

# External object files for target ClippingFunctions
ClippingFunctions_EXTERNAL_OBJECTS =

src/clipping/libClippingFunctions.a: src/clipping/CMakeFiles/ClippingFunctions.dir/ClippingOp.cpp.o
src/clipping/libClippingFunctions.a: src/clipping/CMakeFiles/ClippingFunctions.dir/ReadClipper.cpp.o
src/clipping/libClippingFunctions.a: src/clipping/CMakeFiles/ClippingFunctions.dir/build.make
src/clipping/libClippingFunctions.a: src/clipping/CMakeFiles/ClippingFunctions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libClippingFunctions.a"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && $(CMAKE_COMMAND) -P CMakeFiles/ClippingFunctions.dir/cmake_clean_target.cmake
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ClippingFunctions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/clipping/CMakeFiles/ClippingFunctions.dir/build: src/clipping/libClippingFunctions.a
.PHONY : src/clipping/CMakeFiles/ClippingFunctions.dir/build

src/clipping/CMakeFiles/ClippingFunctions.dir/clean:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping && $(CMAKE_COMMAND) -P CMakeFiles/ClippingFunctions.dir/cmake_clean.cmake
.PHONY : src/clipping/CMakeFiles/ClippingFunctions.dir/clean

src/clipping/CMakeFiles/ClippingFunctions.dir/depend:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bigdreamerxixi/Mutect2Cpp-master /Users/bigdreamerxixi/Mutect2Cpp-master/src/clipping /Users/bigdreamerxixi/Mutect2Cpp-master/build /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/clipping/CMakeFiles/ClippingFunctions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/clipping/CMakeFiles/ClippingFunctions.dir/depend

