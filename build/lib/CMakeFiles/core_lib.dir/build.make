# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

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
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.27.6/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.27.6/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/domsoos/Desktop/research/hipsters/optimization/refactored

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build

# Include any dependencies generated for this target.
include lib/CMakeFiles/core_lib.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/CMakeFiles/core_lib.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/CMakeFiles/core_lib.dir/progress.make

# Include the compile flags for this target's objects.
include lib/CMakeFiles/core_lib.dir/flags.make

lib/CMakeFiles/core_lib.dir/utility.cpp.o: lib/CMakeFiles/core_lib.dir/flags.make
lib/CMakeFiles/core_lib.dir/utility.cpp.o: /Users/domsoos/Desktop/research/hipsters/optimization/refactored/lib/utility.cpp
lib/CMakeFiles/core_lib.dir/utility.cpp.o: lib/CMakeFiles/core_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/CMakeFiles/core_lib.dir/utility.cpp.o"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/CMakeFiles/core_lib.dir/utility.cpp.o -MF CMakeFiles/core_lib.dir/utility.cpp.o.d -o CMakeFiles/core_lib.dir/utility.cpp.o -c /Users/domsoos/Desktop/research/hipsters/optimization/refactored/lib/utility.cpp

lib/CMakeFiles/core_lib.dir/utility.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/core_lib.dir/utility.cpp.i"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/domsoos/Desktop/research/hipsters/optimization/refactored/lib/utility.cpp > CMakeFiles/core_lib.dir/utility.cpp.i

lib/CMakeFiles/core_lib.dir/utility.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/core_lib.dir/utility.cpp.s"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/domsoos/Desktop/research/hipsters/optimization/refactored/lib/utility.cpp -o CMakeFiles/core_lib.dir/utility.cpp.s

lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o: lib/CMakeFiles/core_lib.dir/flags.make
lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o: /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/test_functions.cpp
lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o: lib/CMakeFiles/core_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o -MF CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o.d -o CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o -c /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/test_functions.cpp

lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/core_lib.dir/__/src/test_functions.cpp.i"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/test_functions.cpp > CMakeFiles/core_lib.dir/__/src/test_functions.cpp.i

lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/core_lib.dir/__/src/test_functions.cpp.s"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/test_functions.cpp -o CMakeFiles/core_lib.dir/__/src/test_functions.cpp.s

lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o: lib/CMakeFiles/core_lib.dir/flags.make
lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o: /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/optimization.cpp
lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o: lib/CMakeFiles/core_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o -MF CMakeFiles/core_lib.dir/__/src/optimization.cpp.o.d -o CMakeFiles/core_lib.dir/__/src/optimization.cpp.o -c /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/optimization.cpp

lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/core_lib.dir/__/src/optimization.cpp.i"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/optimization.cpp > CMakeFiles/core_lib.dir/__/src/optimization.cpp.i

lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/core_lib.dir/__/src/optimization.cpp.s"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/optimization.cpp -o CMakeFiles/core_lib.dir/__/src/optimization.cpp.s

lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o: lib/CMakeFiles/core_lib.dir/flags.make
lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o: /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/genetic.cpp
lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o: lib/CMakeFiles/core_lib.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o -MF CMakeFiles/core_lib.dir/__/src/genetic.cpp.o.d -o CMakeFiles/core_lib.dir/__/src/genetic.cpp.o -c /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/genetic.cpp

lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/core_lib.dir/__/src/genetic.cpp.i"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/genetic.cpp > CMakeFiles/core_lib.dir/__/src/genetic.cpp.i

lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/core_lib.dir/__/src/genetic.cpp.s"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/domsoos/Desktop/research/hipsters/optimization/refactored/src/genetic.cpp -o CMakeFiles/core_lib.dir/__/src/genetic.cpp.s

# Object files for target core_lib
core_lib_OBJECTS = \
"CMakeFiles/core_lib.dir/utility.cpp.o" \
"CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o" \
"CMakeFiles/core_lib.dir/__/src/optimization.cpp.o" \
"CMakeFiles/core_lib.dir/__/src/genetic.cpp.o"

# External object files for target core_lib
core_lib_EXTERNAL_OBJECTS =

lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/utility.cpp.o
lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/__/src/test_functions.cpp.o
lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/__/src/optimization.cpp.o
lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/__/src/genetic.cpp.o
lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/build.make
lib/libcore_lib.a: lib/CMakeFiles/core_lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libcore_lib.a"
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/core_lib.dir/cmake_clean_target.cmake
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/core_lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/CMakeFiles/core_lib.dir/build: lib/libcore_lib.a
.PHONY : lib/CMakeFiles/core_lib.dir/build

lib/CMakeFiles/core_lib.dir/clean:
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib && $(CMAKE_COMMAND) -P CMakeFiles/core_lib.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/core_lib.dir/clean

lib/CMakeFiles/core_lib.dir/depend:
	cd /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/domsoos/Desktop/research/hipsters/optimization/refactored /Users/domsoos/Desktop/research/hipsters/optimization/refactored/lib /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib /Users/domsoos/Desktop/research/hipsters/optimization/refactored/build/lib/CMakeFiles/core_lib.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : lib/CMakeFiles/core_lib.dir/depend
