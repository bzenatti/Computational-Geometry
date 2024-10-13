cmake_minimum_required(VERSION 3.1)
cmake_policy(VERSION 3.1)

# Set the project name and version.
project(CGL VERSION 1.0)

# Find required packages.
find_package(CGAL REQUIRED Qt5)

# Build type for performance.
set(CMAKE_BUILD_TYPE Release)

# Set library directories.
set(CGL_SOURCE ${PROJECT_SOURCE_DIR}/src)
set(CGL_INCLUDE ${PROJECT_SOURCE_DIR}/include)
set(CGL_DRAW_SOURCE ${CGL_SOURCE}/Draw)
set(CGL_DRAW_INCLUDE ${CGL_INCLUDE}/Draw)

# Set examples directory.
set(CGL_EXAMPLES ${PROJECT_SOURCE_DIR}/examples)

add_executable(atividade1 ${CGL_EXAMPLES}/atividade1.cpp ${CGL_DRAW_SOURCE}/DrawPoints2.cpp ${CGL_SOURCE}/Random.cpp)
target_include_directories(atividade1 PRIVATE ${CGL_INCLUDE} ${CGL_DRAW_INCLUDE})
target_link_libraries(atividade1 PUBLIC CGAL::CGAL_Basic_viewer)

add_executable(draw_points_2 ${CGL_EXAMPLES}/draw_points_2.cpp ${CGL_DRAW_SOURCE}/DrawPoints2.cpp ${CGL_SOURCE}/Random.cpp)
target_include_directories(draw_points_2 PRIVATE ${CGL_INCLUDE} ${CGL_DRAW_INCLUDE})
target_link_libraries(draw_points_2 PUBLIC CGAL::CGAL_Basic_viewer)

add_executable(draw_vectors_3 ${CGL_EXAMPLES}/draw_vectors_3.cpp ${CGL_DRAW_SOURCE}/DrawVectors3.cpp ${CGL_SOURCE}/Random.cpp)
target_include_directories(draw_vectors_3 PRIVATE ${CGL_INCLUDE} ${CGL_DRAW_INCLUDE})
target_link_libraries(draw_vectors_3 PUBLIC CGAL::CGAL_Basic_viewer)

add_executable(draw_segments_3 ${CGL_EXAMPLES}/draw_segments_3.cpp ${CGL_DRAW_SOURCE}/DrawSegments3.cpp ${CGL_SOURCE}/Random.cpp)
target_include_directories(draw_segments_3 PRIVATE ${CGL_INCLUDE} ${CGL_DRAW_INCLUDE})
target_link_libraries(draw_segments_3 PUBLIC CGAL::CGAL_Basic_viewer)

add_executable(draw_polygon_2 ${CGL_EXAMPLES}/draw_polygon_2.cpp ${CGL_SOURCE}/Random.cpp)
target_include_directories(draw_polygon_2 PRIVATE ${CGL_INCLUDE} ${CGL_DRAW_INCLUDE})
target_link_libraries(draw_polygon_2 PUBLIC CGAL::CGAL_Basic_viewer)

