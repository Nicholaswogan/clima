
add_cython_target(_clima _clima.pyx)

if(WIN32)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMS_WIN64")
endif()

add_library(_clima MODULE ${_clima})
set(_clima_deps
  WaterAdiabatClimate.pyx
)
foreach(PDEP ${_clima_deps})
  set_source_files_properties(_clima.pyx PROPERTIES OBJECT_DEPENDS ${PDEP})
endforeach()

python_extension_module(_clima)

target_link_libraries(_clima clima_wrapper)
target_include_directories(_clima PUBLIC ${NumPy_INCLUDE_DIRS})