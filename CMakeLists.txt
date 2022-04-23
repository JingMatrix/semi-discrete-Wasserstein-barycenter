# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)

project(semi-discrete-Wasserstein-barycenter)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)
set(CMAKE_BUILD_TYPE "Release")
set(CGAL_DATA_DIR ".")

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core Qt5)

include_directories(${CGAL_INCLUDE_DIRS})

set(EXECUTABLE_NAME "semi-discrete")
add_executable(${EXECUTABLE_NAME} "semi-discrete.cpp")
target_link_libraries(
  ${EXECUTABLE_NAME} ${CGAL_LIBRARIES}
  CGAL::CGAL_Basic_viewer
)
