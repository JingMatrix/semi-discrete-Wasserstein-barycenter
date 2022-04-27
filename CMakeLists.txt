# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)

project(semi-discrete-Wasserstein-barycenter)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE "Debug")
set(CGAL_DATA_DIR ".")

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core Qt5)
find_package(GSL REQUIRED)
find_package(GLPK)
include_directories(
  ${CGAL_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${GLPK_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src
)

add_library(power-diagram src/power-diagram.cpp)
add_library(integration src/integration.cpp src/laguerre-cell.cpp)

add_executable(draw-power-diagram src/draw_example.cpp)
add_executable(test src/test.cpp)
add_executable(semi-discrete src/semi-discrete.cpp)

target_link_libraries(
  semi-discrete
  power-diagram
  ${GLPK_LIBRARIES}
)
target_link_libraries(
  integration
  power-diagram ${GSL_LIBRARIES}
)
target_link_libraries(
  draw-power-diagram
  power-diagram
  ${CGAL_LIBRARIES}
  CGAL::CGAL_Basic_viewer
)
target_link_libraries(
  test
  integration
  power-diagram
  ${CGAL_LIBRARIES}
)
