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
CMAKE_SOURCE_DIR = /home/brunopc/Desktop/estudos/Computational-Geometry

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/brunopc/Desktop/estudos/Computational-Geometry/build

# Include any dependencies generated for this target.
include CMakeFiles/polygon_2_triangulate.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/polygon_2_triangulate.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/polygon_2_triangulate.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/polygon_2_triangulate.dir/flags.make

CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o: CMakeFiles/polygon_2_triangulate.dir/flags.make
CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o: ../examples/polygon_2_triangulate.cpp
CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o: CMakeFiles/polygon_2_triangulate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o -MF CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o.d -o CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o -c /home/brunopc/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate.cpp

CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunopc/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate.cpp > CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.i

CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunopc/Desktop/estudos/Computational-Geometry/examples/polygon_2_triangulate.cpp -o CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.s

CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o: CMakeFiles/polygon_2_triangulate.dir/flags.make
CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o: ../src/IO.cpp
CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o: CMakeFiles/polygon_2_triangulate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o -MF CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o.d -o CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o -c /home/brunopc/Desktop/estudos/Computational-Geometry/src/IO.cpp

CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunopc/Desktop/estudos/Computational-Geometry/src/IO.cpp > CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.i

CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunopc/Desktop/estudos/Computational-Geometry/src/IO.cpp -o CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.s

CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o: CMakeFiles/polygon_2_triangulate.dir/flags.make
CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o: ../src/Partition.cpp
CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o: CMakeFiles/polygon_2_triangulate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o -MF CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o.d -o CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o -c /home/brunopc/Desktop/estudos/Computational-Geometry/src/Partition.cpp

CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunopc/Desktop/estudos/Computational-Geometry/src/Partition.cpp > CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.i

CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunopc/Desktop/estudos/Computational-Geometry/src/Partition.cpp -o CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.s

CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o: CMakeFiles/polygon_2_triangulate.dir/flags.make
CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o: ../src/Primitives.cpp
CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o: CMakeFiles/polygon_2_triangulate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o -MF CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o.d -o CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o -c /home/brunopc/Desktop/estudos/Computational-Geometry/src/Primitives.cpp

CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunopc/Desktop/estudos/Computational-Geometry/src/Primitives.cpp > CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.i

CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunopc/Desktop/estudos/Computational-Geometry/src/Primitives.cpp -o CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.s

CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o: CMakeFiles/polygon_2_triangulate.dir/flags.make
CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o: ../src/Draw/DrawSegments3.cpp
CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o: CMakeFiles/polygon_2_triangulate.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o -MF CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o.d -o CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o -c /home/brunopc/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp

CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/brunopc/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp > CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.i

CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/brunopc/Desktop/estudos/Computational-Geometry/src/Draw/DrawSegments3.cpp -o CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.s

# Object files for target polygon_2_triangulate
polygon_2_triangulate_OBJECTS = \
"CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o" \
"CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o" \
"CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o" \
"CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o" \
"CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o"

# External object files for target polygon_2_triangulate
polygon_2_triangulate_EXTERNAL_OBJECTS =

polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/examples/polygon_2_triangulate.cpp.o
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/src/IO.cpp.o
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/src/Partition.cpp.o
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/src/Primitives.cpp.o
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/src/Draw/DrawSegments3.cpp.o
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/build.make
polygon_2_triangulate: libCGAL_Qt5_moc_and_resources.a
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libgmpxx.so
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libmpfr.so
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libgmp.so
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libQt5Svg.so.5.15.3
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.15.3
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.3
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.3
polygon_2_triangulate: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.3
polygon_2_triangulate: CMakeFiles/polygon_2_triangulate.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable polygon_2_triangulate"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polygon_2_triangulate.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/polygon_2_triangulate.dir/build: polygon_2_triangulate
.PHONY : CMakeFiles/polygon_2_triangulate.dir/build

CMakeFiles/polygon_2_triangulate.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/polygon_2_triangulate.dir/cmake_clean.cmake
.PHONY : CMakeFiles/polygon_2_triangulate.dir/clean

CMakeFiles/polygon_2_triangulate.dir/depend:
	cd /home/brunopc/Desktop/estudos/Computational-Geometry/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/brunopc/Desktop/estudos/Computational-Geometry /home/brunopc/Desktop/estudos/Computational-Geometry /home/brunopc/Desktop/estudos/Computational-Geometry/build /home/brunopc/Desktop/estudos/Computational-Geometry/build /home/brunopc/Desktop/estudos/Computational-Geometry/build/CMakeFiles/polygon_2_triangulate.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/polygon_2_triangulate.dir/depend

