cmake_minimum_required(VERSION "3.14")

project(CLIMA LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
endif()

# includes
include(cmake/CPM.cmake)

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

