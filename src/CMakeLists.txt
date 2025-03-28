add_library(differentia
  differentia_const.f90
  differentia_dual.f90
  differentia_derivative.f90
  differentia.f90
)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(differentia PRIVATE -Wunused -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(differentia PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()