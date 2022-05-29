
set(CLIMA_TESTS test_clima test_adiabat test_kasting)

foreach(CLIMA_TEST ${CLIMA_TESTS})
  add_executable(${CLIMA_TEST} "${CLIMA_TEST}.f90")
  target_link_libraries(${CLIMA_TEST} clima)
  target_include_directories(${CLIMA_TEST} PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})
endforeach()
