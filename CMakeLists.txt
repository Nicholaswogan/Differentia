cmake_minimum_required(VERSION "3.14" )

project(DIFFERENTIA LANGUAGES Fortran VERSION "0.1.4")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

add_subdirectory(src)
add_subdirectory(test)

