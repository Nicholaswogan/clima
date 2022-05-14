
add_executable(test_clima test_clima.f90)
target_link_libraries(test_clima clima)
target_include_directories(test_clima PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(validate validate.f90)
target_link_libraries(validate clima)
target_include_directories(validate PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})