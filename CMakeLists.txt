cmake_minimum_required(VERSION "3.14")

project(CLIMA LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# includes
include(cmake/CPM.cmake)

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

