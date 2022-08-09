
add_library(clima_wrapper 
  wrapper_utils.f90 
  WaterAdiabatClimate_wrapper.f90 
)
target_link_libraries(clima_wrapper clima)
