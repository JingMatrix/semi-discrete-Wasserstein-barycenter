# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)

project(semi-discrete-Wasserstein-barycenter)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE "Debug")
set(CGAL_DATA_DIR ".")

option(USE_EXACT_KERNEL "whether use the CGAL exact kernel" OFF)
if(USE_EXACT_KERNEL)
  add_definitions(-DUSE_EXACT_KERNEL)
endif(USE_EXACT_KERNEL)

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core Qt5)
find_package(GSL REQUIRED)
find_package(GLPK)
include_directories(
  ${CGAL_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${GLPK_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src
)

add_library(
  power-diagram
  src/power-diagram.cpp
  src/integration.cpp
  src/crop-cell.cpp
  src/crop.cpp
  src/draw.cpp
)
add_library(
  barycenter
  src/marginals.cpp
  src/numerical_solver.cpp
  src/linear-programming.cpp
)

add_executable(draw-power-diagram src/qt_draw_example.cpp)
add_executable(test src/test.cpp)

target_link_libraries(
  power-diagram
  ${CGAL_LIBRARIES}
)
target_link_libraries(
  barycenter
  ${GLPK_LIBRARIES}
)
target_link_libraries(
  draw-power-diagram
  power-diagram
  CGAL::CGAL_Basic_viewer
)
target_link_libraries(
  test
  barycenter
  power-diagram
)
