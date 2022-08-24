cmake_minimum_required(VERSION "3.14")

project(CLIMA LANGUAGES Fortran C VERSION "0.2.1")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# includes
include(cmake/CPM.cmake)

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

option(BUILD_PYTHON_CLIMA "if ON, then will build a python
version via Cython" OFF)

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

if (BUILD_PYTHON_CLIMA)
  if (NOT SKBUILD)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  
  add_subdirectory(clima)
  
  if (SKBUILD)
    install(TARGETS _clima DESTINATION clima)
  else()
    install(TARGETS _clima DESTINATION ${CMAKE_SOURCE_DIR}/clima)
  endif()
                
endif()

