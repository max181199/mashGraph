# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/max/work/Graph/task1_src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/max/work/Graph/task1_src/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/rt.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/rt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rt.dir/flags.make

CMakeFiles/rt.dir/main.cpp.o: CMakeFiles/rt.dir/flags.make
CMakeFiles/rt.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/max/work/Graph/task1_src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rt.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rt.dir/main.cpp.o -c /Users/max/work/Graph/task1_src/main.cpp

CMakeFiles/rt.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rt.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/max/work/Graph/task1_src/main.cpp > CMakeFiles/rt.dir/main.cpp.i

CMakeFiles/rt.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rt.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/max/work/Graph/task1_src/main.cpp -o CMakeFiles/rt.dir/main.cpp.s

CMakeFiles/rt.dir/Bitmap.cpp.o: CMakeFiles/rt.dir/flags.make
CMakeFiles/rt.dir/Bitmap.cpp.o: ../Bitmap.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/max/work/Graph/task1_src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/rt.dir/Bitmap.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/rt.dir/Bitmap.cpp.o -c /Users/max/work/Graph/task1_src/Bitmap.cpp

CMakeFiles/rt.dir/Bitmap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rt.dir/Bitmap.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/max/work/Graph/task1_src/Bitmap.cpp > CMakeFiles/rt.dir/Bitmap.cpp.i

CMakeFiles/rt.dir/Bitmap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rt.dir/Bitmap.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/max/work/Graph/task1_src/Bitmap.cpp -o CMakeFiles/rt.dir/Bitmap.cpp.s

# Object files for target rt
rt_OBJECTS = \
"CMakeFiles/rt.dir/main.cpp.o" \
"CMakeFiles/rt.dir/Bitmap.cpp.o"

# External object files for target rt
rt_EXTERNAL_OBJECTS =

rt: CMakeFiles/rt.dir/main.cpp.o
rt: CMakeFiles/rt.dir/Bitmap.cpp.o
rt: CMakeFiles/rt.dir/build.make
rt: CMakeFiles/rt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/max/work/Graph/task1_src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable rt"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/rt.dir/build: rt

.PHONY : CMakeFiles/rt.dir/build

CMakeFiles/rt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rt.dir/clean

CMakeFiles/rt.dir/depend:
	cd /Users/max/work/Graph/task1_src/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/max/work/Graph/task1_src /Users/max/work/Graph/task1_src /Users/max/work/Graph/task1_src/cmake-build-debug /Users/max/work/Graph/task1_src/cmake-build-debug /Users/max/work/Graph/task1_src/cmake-build-debug/CMakeFiles/rt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rt.dir/depend

