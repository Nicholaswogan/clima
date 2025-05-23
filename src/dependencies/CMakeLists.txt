
CPMAddPackage(
  NAME differentia
  VERSION 0.1.4
  OPTIONS 
    "BUILD_EXECUTABLE OFF"
  GITHUB_REPOSITORY "Nicholaswogan/Differentia"
  GIT_TAG "v0.1.4"
  EXCLUDE_FROM_ALL ON
)

# fortran-yaml-c
CPMAddPackage(
  NAME fortran-yaml-c
  VERSION 0.2.6
  GITHUB_REPOSITORY "Nicholaswogan/fortran-yaml-c"
  GIT_TAG "v0.2.6"
  EXCLUDE_FROM_ALL ON
)

# h5fortran
# Here, I've made a tiny modification to v4.5.0 in a fork of the
# original library (geospace-code/h5fortran) so that it does not 
# exclude conda libraries. 
CPMAddPackage(
  NAME h5fortran
  VERSION 4.5.0_c55f15e94300be746b22b4b840a6c4051d9fb323
  OPTIONS
    "CMAKE_POSITION_INDEPENDENT_CODE ON"
  GITHUB_REPOSITORY "Nicholaswogan/h5fortran"
  GIT_TAG "c55f15e94300be746b22b4b840a6c4051d9fb323"
  EXCLUDE_FROM_ALL ON
)

# futils
CPMAddPackage(
  NAME futils
  VERSION 0.1.14
  GITHUB_REPOSITORY "Nicholaswogan/futils"
  GIT_TAG "v0.1.14"
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
target_link_libraries(finterp differentia)

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
