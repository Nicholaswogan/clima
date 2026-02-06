# First download data, if needed
if (NOT EXISTS "${CMAKE_SOURCE_DIR}/data")
  file(DOWNLOAD 
    "https://github.com/Nicholaswogan/photochem_clima_data/archive/v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
    "photochem_clima_data_v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
  )
  file(ARCHIVE_EXTRACT 
    INPUT "photochem_clima_data_v${PHOTOCHEM_CLIMA_DATA_VERSION}.tar.gz"
    DESTINATION ${CMAKE_BINARY_DIR}
  )
  file(RENAME 
    "${CMAKE_BINARY_DIR}/photochem_clima_data-${PHOTOCHEM_CLIMA_DATA_VERSION}/photochem_clima_data/data"
    "${CMAKE_SOURCE_DIR}/data"
  )
endif()

set(CLIMA_TESTS test_radtran test_adiabat)

if (BUILD_EXECUTABLES)
  foreach(CLIMA_TEST ${CLIMA_TESTS})
    add_executable(${CLIMA_TEST} "${CLIMA_TEST}.f90")
    target_link_libraries(${CLIMA_TEST} clima)
    target_include_directories(${CLIMA_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
  endforeach()
endif()