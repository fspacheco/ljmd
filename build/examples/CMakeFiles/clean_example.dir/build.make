# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/massimo/MHPC/ljmd

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/massimo/MHPC/ljmd/build

# Utility rule file for clean_example.

# Include any custom commands dependencies for this target.
include examples/CMakeFiles/clean_example.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/clean_example.dir/progress.make

examples/CMakeFiles/clean_example:
	cd /home/massimo/MHPC/ljmd/build/examples && rm -f *.xyz *.dat

clean_example: examples/CMakeFiles/clean_example
clean_example: examples/CMakeFiles/clean_example.dir/build.make
.PHONY : clean_example

# Rule to build all files generated by this target.
examples/CMakeFiles/clean_example.dir/build: clean_example
.PHONY : examples/CMakeFiles/clean_example.dir/build

examples/CMakeFiles/clean_example.dir/clean:
	cd /home/massimo/MHPC/ljmd/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/clean_example.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/clean_example.dir/clean

examples/CMakeFiles/clean_example.dir/depend:
	cd /home/massimo/MHPC/ljmd/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/massimo/MHPC/ljmd /home/massimo/MHPC/ljmd/examples /home/massimo/MHPC/ljmd/build /home/massimo/MHPC/ljmd/build/examples /home/massimo/MHPC/ljmd/build/examples/CMakeFiles/clean_example.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/clean_example.dir/depend

