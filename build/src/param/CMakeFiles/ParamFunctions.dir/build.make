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
include src/param/CMakeFiles/ParamFunctions.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/param/CMakeFiles/ParamFunctions.dir/compiler_depend.make

# Include the progress variables for this target.
include src/param/CMakeFiles/ParamFunctions.dir/progress.make

# Include the compile flags for this target's objects.
include src/param/CMakeFiles/ParamFunctions.dir/flags.make

src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o: src/param/CMakeFiles/ParamFunctions.dir/flags.make
src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o: ../src/param/ParamUtils.cpp
src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o: src/param/CMakeFiles/ParamFunctions.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o -MF CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o.d -o CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o -c /Users/bigdreamerxixi/Mutect2Cpp-master/src/param/ParamUtils.cpp

src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.i"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/bigdreamerxixi/Mutect2Cpp-master/src/param/ParamUtils.cpp > CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.i

src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.s"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/bigdreamerxixi/Mutect2Cpp-master/src/param/ParamUtils.cpp -o CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.s

# Object files for target ParamFunctions
ParamFunctions_OBJECTS = \
"CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o"

# External object files for target ParamFunctions
ParamFunctions_EXTERNAL_OBJECTS =

src/param/libParamFunctions.a: src/param/CMakeFiles/ParamFunctions.dir/ParamUtils.cpp.o
src/param/libParamFunctions.a: src/param/CMakeFiles/ParamFunctions.dir/build.make
src/param/libParamFunctions.a: src/param/CMakeFiles/ParamFunctions.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/bigdreamerxixi/Mutect2Cpp-master/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libParamFunctions.a"
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && $(CMAKE_COMMAND) -P CMakeFiles/ParamFunctions.dir/cmake_clean_target.cmake
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ParamFunctions.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/param/CMakeFiles/ParamFunctions.dir/build: src/param/libParamFunctions.a
.PHONY : src/param/CMakeFiles/ParamFunctions.dir/build

src/param/CMakeFiles/ParamFunctions.dir/clean:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param && $(CMAKE_COMMAND) -P CMakeFiles/ParamFunctions.dir/cmake_clean.cmake
.PHONY : src/param/CMakeFiles/ParamFunctions.dir/clean

src/param/CMakeFiles/ParamFunctions.dir/depend:
	cd /Users/bigdreamerxixi/Mutect2Cpp-master/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bigdreamerxixi/Mutect2Cpp-master /Users/bigdreamerxixi/Mutect2Cpp-master/src/param /Users/bigdreamerxixi/Mutect2Cpp-master/build /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param /Users/bigdreamerxixi/Mutect2Cpp-master/build/src/param/CMakeFiles/ParamFunctions.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/param/CMakeFiles/ParamFunctions.dir/depend

