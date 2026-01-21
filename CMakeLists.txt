cmake_minimum_required(VERSION "3.14")
cmake_policy(SET CMP0148 OLD) # To suppress a warning emerging from scikit-build

project(Clima LANGUAGES Fortran C VERSION "0.5.15")
set(PHOTOCHEM_CLIMA_DATA_VERSION "0.3.1")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# includes
include(cmake/CPM.cmake)
include(cmake/fypp.cmake)

# options
option(SKBUILD "Should be ON of being build by skbuild, 
and OFF of being build by regular cmake" OFF)
message(STATUS "The project is built using scikit-build: ${SKBUILD}")

option(BUILD_EXECUTABLES "if ON, then will build the
Fortran executables" ON)

option(BUILD_PYTHON_CLIMA "if ON, then will build a python
version via Cython" OFF)

if (NOT DEFINED PYTHON_CLIMA_DESTINATION)
  set(PYTHON_CLIMA_DESTINATION "clima" CACHE STRING 
  "folder that python version of climate is installed")
endif()

# src
add_subdirectory(src)

# test
add_subdirectory(tests)

if (BUILD_PYTHON_CLIMA)
  if (NOT SKBUILD)
    if (NOT DEFINED SKBUILD_CMAKE_MODULE_DIR)
      # Here, we try to find scikit-build cmake modules
      find_package(Python COMPONENTS Development)
      set(SKBUILD_CMAKE_MODULE_DIR "${Python_LIBRARY_DIRS}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/skbuild/resources/cmake")
    endif()
    if (EXISTS ${SKBUILD_CMAKE_MODULE_DIR})
      message(STATUS "Scikit-build CMake modules: ${SKBUILD_CMAKE_MODULE_DIR}")
    else()
      message(FATAL_ERROR "Failed to find scikit-build CMake modules in directory: ${SKBUILD_CMAKE_MODULE_DIR}")
    endif()
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${SKBUILD_CMAKE_MODULE_DIR})
  endif()
  
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)
  find_package(Cython REQUIRED)
  
  add_subdirectory(clima)
  
  if (SKBUILD)
    install(TARGETS _clima DESTINATION ${PYTHON_CLIMA_DESTINATION})
  else()
    install(TARGETS _clima DESTINATION ${CMAKE_SOURCE_DIR}/${PYTHON_CLIMA_DESTINATION})
  endif()
                
endif()

