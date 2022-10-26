
add_library(clima_c_api 
  clima_c_api.f90 
)
set(clima_c_api_deps
  futils.f90
  OpticalProperties.f90
  ClimaRadtranWrk.f90
  Radtran.f90
  WaterAdiabatClimate.f90
)
foreach(PDEP ${_clima_deps})
  set_source_files_properties(clima_c_api.f90 PROPERTIES OBJECT_DEPENDS ${PDEP})
endforeach()
target_link_libraries(clima_c_api clima)
