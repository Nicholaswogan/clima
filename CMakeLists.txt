cmake_minimum_required(VERSION "3.14")

project(CLIMA LANGUAGES Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wunused")
  
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finline-functions -ftree-vectorize -funroll-loops -march=native")
  endif()
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

# includes
include(cmake/CPM.cmake)

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

