# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/brunoz/Desktop/estudos/Computational-Geometry

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/brunoz/Desktop/estudos/Computational-Geometry/build

# Include any dependencies generated for this target.
include CMakeFiles/polygon_2_triangulate_good.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/polygon_2_triangulate_good.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/polygon_2_triangulate_good.dir/flags.make

CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o: CMakeFiles/polygon_2_triangulate_good.dir/flags.make
CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o: ../examples/polygon_2_triangulate_good.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o -c /home/brunoz/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate_good.cpp

CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunoz/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate_good.cpp > CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.i

CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunoz/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate_good.cpp -o CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.s

CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o: CMakeFiles/polygon_2_triangulate_good.dir/flags.make
CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o: ../src/IO.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o -c /home/brunoz/Desktop/estudos/Computational-Geometry/src/IO.cpp

CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunoz/Desktop/estudos/Computational-Geometry/src/IO.cpp > CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.i

CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunoz/Desktop/estudos/Computational-Geometry/src/IO.cpp -o CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.s

CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o: CMakeFiles/polygon_2_triangulate_good.dir/flags.make
CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o: ../src/Partition.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o -c /home/brunoz/Desktop/estudos/Computational-Geometry/src/Partition.cpp

CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunoz/Desktop/estudos/Computational-Geometry/src/Partition.cpp > CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.i

CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunoz/Desktop/estudos/Computational-Geometry/src/Partition.cpp -o CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.s

CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o: CMakeFiles/polygon_2_triangulate_good.dir/flags.make
CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o: ../src/Primitives.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o -c /home/brunoz/Desktop/estudos/Computational-Geometry/src/Primitives.cpp

CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunoz/Desktop/estudos/Computational-Geometry/src/Primitives.cpp > CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.i

CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunoz/Desktop/estudos/Computational-Geometry/src/Primitives.cpp -o CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.s

CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o: CMakeFiles/polygon_2_triangulate_good.dir/flags.make
CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o: ../src/Draw/DrawSegments3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o -c /home/brunoz/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp

CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunoz/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp > CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.i

CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunoz/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp -o CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.s

# Object files for target polygon_2_triangulate_good
polygon_2_triangulate_good_OBJECTS = \
"CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o" \
"CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o" \
"CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o" \
"CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o" \
"CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o"

# External object files for target polygon_2_triangulate_good
polygon_2_triangulate_good_EXTERNAL_OBJECTS =

polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/examples/polygon_2_triangulate_good.cpp.o
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/src/IO.cpp.o
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/src/Partition.cpp.o
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/src/Primitives.cpp.o
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/src/Draw/DrawSegments3.cpp.o
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/build.make
polygon_2_triangulate_good: libCGAL_Qt5_moc_and_resources.a
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libgmpxx.so
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libmpfr.so
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libgmp.so
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libQt5Svg.so.5.12.8
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.12.8
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.12.8
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.12.8
polygon_2_triangulate_good: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.12.8
polygon_2_triangulate_good: CMakeFiles/polygon_2_triangulate_good.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable polygon_2_triangulate_good"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polygon_2_triangulate_good.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/polygon_2_triangulate_good.dir/build: polygon_2_triangulate_good

.PHONY : CMakeFiles/polygon_2_triangulate_good.dir/build

CMakeFiles/polygon_2_triangulate_good.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/polygon_2_triangulate_good.dir/cmake_clean.cmake
.PHONY : CMakeFiles/polygon_2_triangulate_good.dir/clean

CMakeFiles/polygon_2_triangulate_good.dir/depend:
	cd /home/brunoz/Desktop/estudos/Computational-Geometry/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/brunoz/Desktop/estudos/Computational-Geometry /home/brunoz/Desktop/estudos/Computational-Geometry /home/brunoz/Desktop/estudos/Computational-Geometry/build /home/brunoz/Desktop/estudos/Computational-Geometry/build /home/brunoz/Desktop/estudos/Computational-Geometry/build/CMakeFiles/polygon_2_triangulate_good.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/polygon_2_triangulate_good.dir/depend

