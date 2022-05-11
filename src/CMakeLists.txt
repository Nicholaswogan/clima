
add_subdirectory(dependencies)

set(CLIMATE_SOURCES
  climate/clima_climate_rhs.f90
  climate/clima_climate_integrate.f90
  climate/clima_climate.f90
)

add_library(clima 
  clima_const.f90
  clima_types.f90
  clima_eqns.f90
  clima_input.f90
  clima_twostream.f90
  clima_radtran.f90
  ${CLIMATE_SOURCES}
  clima.f90
)
target_link_libraries(clima fortran-yaml-c finterp h5fortran futils fortran_stdlib mrgrnk sundials_fcvode_mod sundials_cvode dop853)

# Add special compiler flags to Clima
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(clima PRIVATE -Wunused)
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(clima PRIVATE -finline-functions -ftree-vectorize -funroll-loops -march=native)
  endif()
endif()

option(BUILD_WITH_OPENMP "Compile with muli-threading" OFF)
if (BUILD_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_Fortran_FOUND)
    target_compile_options(clima PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(clima PUBLIC ${OpenMP_Fortran_FLAGS})
  endif()
endif()
message (STATUS "Building with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

