
add_library(clima 
  clima_const.f90
  clima_types.f90
  clima.f90
)
target_link_libraries(clima fortran-yaml-c)