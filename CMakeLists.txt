cmake_minimum_required(VERSION "3.14")

project(CLIMA LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# dependencies
include(cmake/CPM.cmake)
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.1.2
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.1.2"
  EXCLUDE_FROM_ALL ON
)

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

