
add_subdirectory(dependencies)

set(BASE_SOURCES
  clima_const.f90
  clima_eqns.f90
  clima_types.f90
  clima_types_create.f90
)

set(RADIATE_SOURCES
  radtran/clima_radtran_types.f90
  radtran/clima_radtran_types_create.f90
  radtran/clima_radtran_twostream.f90
  radtran/clima_radtran_radiate.f90
  radtran/clima_radtran.f90
)

set(CLIMATE_SOURCES
  climate/clima_climate.f90
  climate/clima_climate_rhs.f90
  climate/clima_climate_integrate.f90
)

add_library(clima 
  ${BASE_SOURCES}
  ${RADIATE_SOURCES}
  ${CLIMATE_SOURCES}
)
target_link_libraries(clima fortran-yaml-c finterp h5fortran futils fortran_stdlib mrgrnk dop853 minpack)

# Add special compiler flags to Clima
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(clima PRIVATE -Wunused -fcheck=all,no-array-temps -fimplicit-none)
  endif()
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(clima PRIVATE -finline-functions -ftree-vectorize -funroll-loops)
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

