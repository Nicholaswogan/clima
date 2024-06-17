
add_subdirectory(dependencies)
# Find lapack
find_package(LAPACK REQUIRED)

configure_file(clima_version.f90.in ${CMAKE_CURRENT_BINARY_DIR}/clima_version.f90)

set(fppFiles 
  clima_saturationdata.f90
)

fypp_f90("" "${fppFiles}" outFiles)

set(BASE_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/clima_version.f90
  clima_const.f90
  clima_useful.f90
  clima_eqns.f90
  clima_eqns_water.f90
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

set(ADIABAT_SOURCES
  adiabat/clima_adiabat_general.f90
  adiabat/clima_adiabat_rc.f90
  adiabat/clima_adiabat_solve.f90
  adiabat/clima_adiabat.f90
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
  ${ADIABAT_SOURCES}
  ${outFiles}
  clima.f90
)
target_link_libraries(clima 
  fortran-yaml-c  
  h5fortran 
  futils 
  finterp
  dop853
  minpack
  forwarddiff
  ${LAPACK_LIBRARIES}
)

# This is important for making Python wrapper
set_property(TARGET clima PROPERTY POSITION_INDEPENDENT_CODE ON)

# Add special compiler flags to Clima
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(clima PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(clima PRIVATE -fcheck=all,no-array-temps)
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
message (STATUS "Building clima with OpenMP multi-threading = ${BUILD_WITH_OPENMP}")

