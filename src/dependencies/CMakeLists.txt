
CPMAddPackage(
  NAME forwarddiff
  VERSION 0.1.2
  OPTIONS 
    "BUILD_EXECUTABLE OFF"
  GITHUB_REPOSITORY "Nicholaswogan/forwarddiff"
  GIT_TAG "v0.1.2"
  EXCLUDE_FROM_ALL ON
)

# fortran-yaml-c
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.4
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.4"
  EXCLUDE_FROM_ALL ON
)

# h5fortran
CPMAddPackage(
  NAME h5fortran
  VERSION 4.5.0
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
  GITHUB_REPOSITORY "geospace-code/h5fortran"
  GIT_TAG "v4.5.0"
  EXCLUDE_FROM_ALL ON
)

# futils
CPMAddPackage(
  NAME futils
  VERSION 0.1.13
  GITHUB_REPOSITORY "Nicholaswogan/futils"
  GIT_TAG "v0.1.13"
  EXCLUDE_FROM_ALL ON
)

# finterp
# CPMAddPackage(
#   NAME finterp
#   VERSION 1.3.0
#   GITHUB_REPOSITORY "jacobwilliams/finterp"
#   GIT_TAG "1.3.0"
#   DOWNLOAD_ONLY ON
# )
set(fppFiles 
  linear_interpolation_module.F90
)
fypp_f90("" "${fppFiles}" outFiles)
add_library(finterp ${outFiles})
target_link_libraries(finterp forwarddiff)

CPMAddPackage(
  NAME dop853
  VERSION 1.2.0
  GITHUB_REPOSITORY "jacobwilliams/dop853"
  GIT_TAG "1.2.0"
  DOWNLOAD_ONLY ON
)

add_library(dop853 
  ${dop853_SOURCE_DIR}/src/dop853_constants.f90
  ${dop853_SOURCE_DIR}/src/dop853_module.f90
)

CPMAddPackage(
  NAME minpack
  VERSION 2.0.0
  GITHUB_REPOSITORY "fortran-lang/minpack"
  GIT_TAG "v2.0.0-rc.1"
  DOWNLOAD_ONLY ON
)

add_library(minpack
  ${minpack_SOURCE_DIR}/src/minpack.f90
)
