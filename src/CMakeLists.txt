
add_subdirectory(dependencies)

add_library(clima 
  clima_const.f90
  clima_types.f90
  clima_input.f90
  clima_climate.f90
  clima.f90
)
target_link_libraries(clima fortran-yaml-c finterp h5fortran futils)