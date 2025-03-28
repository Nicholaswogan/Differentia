add_executable(test_differentia test_differentia.f90)
target_link_libraries(test_differentia differentia)
target_include_directories(test_differentia PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(test_sparse test_sparse.f90)
target_link_libraries(test_sparse differentia)
target_include_directories(test_sparse PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

find_program(FYPP fypp REQUIRED)
include(${PROJECT_SOURCE_DIR}/cmake/fypp.cmake)

set(fppFiles fypp_example.fypp)
fypp_f90("" "${fppFiles}" outFiles)
add_executable(fypp_example ${outFiles})
target_link_libraries(fypp_example differentia)
target_include_directories(fypp_example PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})